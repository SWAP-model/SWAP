!> SWAP ensemble: N independent column states sharing read-only config(s),
!! stepped sequentially. Backing store for the XMI coupling kernel. One
!! ensemble per process. Per-column config indirection (column_config ->
!! configs) is built in; the demo uses a single shared config (M=1).
module swap_ensemble_mod
   use iso_fortran_env, only: real64, int32
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_mod,        only: swap_init_from_loaded_config, swap_run_step, swap_close
   use load_swap_config_mod, only: load_swap_config
   use error_mod,       only: error_collection_t, &
                              set_library_mode, library_fatal_raised, clear_library_fatal
   implicit none
   private

   public :: ensemble_init, ensemble_step_day, ensemble_finalize
   public :: ensemble_ncol, ensemble_all_swbotb1
   public :: ensemble_set_gwl, ensemble_qbot_volume, ensemble_storage_coef
   public :: ensemble_current_t1900, ensemble_start_t1900, ensemble_end_t1900
   public :: gwl, qbot_volume, storage_coef

   integer, parameter :: STORAGE_COEF_DEFAULT_X100 = 15   ! sy = 0.15 (smoke placeholder)

   type(swap_state_t),  allocatable, save :: columns(:)
   type(swap_config_t), allocatable, save, target :: configs(:)
   integer,             allocatable, save :: column_config(:)
   integer,             save :: ncol = 0

   real(real64), allocatable, save, target :: gwl(:)          ! MODFLOW head, metres (driver writes)
   real(real64), allocatable, save, target :: qbot_volume(:)  ! recharge depth over step, metres (we write)
   real(real64), allocatable, save, target :: storage_coef(:) ! specific yield, - (we write)

contains

   integer function ensemble_init(config_file, ncol_in) result(rc)
      character(len=*), intent(in) :: config_file
      integer,          intent(in) :: ncol_in
      type(error_collection_t) :: errors
      integer :: i
      rc = 0
      call set_library_mode(.true.)
      call clear_library_fatal()
      ncol = ncol_in

      ! Load + validate + finalize the shared config, mirroring swap_init's
      ! TOML pipeline exactly so the ensemble config is constructed identically
      ! to the standalone path. In library mode a fatal sets a poll-able flag
      ! (abort_if_fatal does not error stop) which we check after.
      allocate(configs(1))
      call load_swap_config(config_file, configs(1), errors)
      call configs(1)%validate(errors)
      call configs(1)%finalize(errors)
      call errors%abort_if_fatal()
      if (library_fatal_raised()) then; rc = 1; return; end if

      allocate(column_config(ncol)); column_config = 1
      allocate(columns(ncol))
      do i = 1, ncol
         call swap_init_from_loaded_config(columns(i), configs(column_config(i)))
         if (library_fatal_raised()) then; rc = 2; return; end if
         columns(i)%soilwater%flcoupled_gwl = .true.
      end do

      allocate(gwl(ncol),          source=0.0_real64)
      allocate(qbot_volume(ncol),  source=0.0_real64)
      allocate(storage_coef(ncol), source=real(STORAGE_COEF_DEFAULT_X100, real64)/100.0_real64)

      if (.not. ensemble_all_swbotb1()) rc = 3
   end function ensemble_init

   integer function ensemble_step_day() result(rc)
      integer      :: i
      real(real64) :: dt_days
      rc = 0
      do i = 1, ncol
         columns(i)%soilwater%gwl_injected = gwl(i) * 100.0_real64   ! m -> cm
         call swap_run_step(columns(i), configs(column_config(i)))
         if (library_fatal_raised()) then; rc = 1; return; end if
         dt_days = columns(i)%timecontrol%dt
         qbot_volume(i) = (-columns(i)%soilwater%qbot) * 0.01_real64 * dt_days   ! cm/d -> m depth, sign flip
      end do
   end function ensemble_step_day

   integer function ensemble_finalize() result(rc)
      integer :: i
      rc = 0
      do i = 1, ncol
         call swap_close(columns(i), configs(column_config(i)))
      end do
      if (allocated(columns))       deallocate(columns)
      if (allocated(configs))       deallocate(configs)
      if (allocated(column_config)) deallocate(column_config)
      if (allocated(gwl))           deallocate(gwl)
      if (allocated(qbot_volume))   deallocate(qbot_volume)
      if (allocated(storage_coef))  deallocate(storage_coef)
      ncol = 0
   end function ensemble_finalize

   integer function ensemble_ncol(); ensemble_ncol = ncol; end function

   logical function ensemble_all_swbotb1()
      integer :: i
      ensemble_all_swbotb1 = (ncol > 0)
      do i = 1, ncol
         if (columns(i)%soilwater%swbotb_runtime /= 1) ensemble_all_swbotb1 = .false.
      end do
   end function

   subroutine ensemble_set_gwl(i, head_m)
      integer,      intent(in) :: i
      real(real64), intent(in) :: head_m
      gwl(i) = head_m
   end subroutine

   real(real64) function ensemble_qbot_volume(i); integer, intent(in) :: i
      ensemble_qbot_volume = qbot_volume(i); end function
   real(real64) function ensemble_storage_coef(i); integer, intent(in) :: i
      ensemble_storage_coef = storage_coef(i); end function

   real(real64) function ensemble_current_t1900()
      ensemble_current_t1900 = columns(1)%timecontrol%t1900
   end function
   real(real64) function ensemble_start_t1900()
      ensemble_start_t1900 = columns(1)%timecontrol%tstart
   end function
   real(real64) function ensemble_end_t1900()
      ensemble_end_t1900 = columns(1)%timecontrol%tend
   end function

end module swap_ensemble_mod
