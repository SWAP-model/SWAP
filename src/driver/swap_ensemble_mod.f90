!> SWAP ensemble: N independent column states sharing read-only config(s),
!! stepped sequentially. Backing store for the XMI coupling kernel. One
!! ensemble per process. Per-column config indirection (column_config ->
!! configs) is built in; the demo uses a single shared config (M=1).
module swap_ensemble_mod
   use iso_fortran_env, only: real64, int32
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_mod,        only: swap_init, swap_init_from_loaded_config, swap_run_step, swap_close
   use soilhydraulics_mod, only: reequilibrate_column_to_gwl
   use load_swap_config_mod, only: load_swap_config
   use error_mod,       only: error_collection_t, &
                              set_library_mode, library_fatal_raised, clear_library_fatal
   use diagnostics_mod, only: diag_overrides_t, default_embedded_config, &
                              read_logging_overrides_from_file, init_logging
   implicit none
   private

   public :: ensemble_init, ensemble_step_day, ensemble_finalize
   public :: ensemble_ncol, ensemble_all_swbotb1
   public :: ensemble_set_gwl, ensemble_qbot_volume, ensemble_storage_coef
   public :: ensemble_current_t1900, ensemble_start_t1900, ensemble_end_t1900
   public :: gwl, qbot_volume, storage_coef
   ! T1-G′ sub-arc 2: single-column standalone mode + shared facade plumbing.
   public :: ensemble_allocate_single, ensemble_init_single, ensemble_is_coupled
   public :: ensemble_column1, ensemble_config1
   public :: read_ncol_sidecar, ensemble_last_error

   integer, parameter :: STORAGE_COEF_DEFAULT_X100 = 15   ! sy = 0.15 (smoke placeholder)

   type(swap_state_t),  allocatable, save, target :: columns(:)
   type(swap_config_t), allocatable, save, target :: configs(:)
   integer,             allocatable, save :: column_config(:)
   integer,             save :: ncol = 0
   logical,             save :: first_step_done = .false.  ! gate per-cell gwl re-seat
   logical,             save :: coupled_mode    = .false.  ! set by ensemble_init (coupled) / cleared by single-mode init
   character(len=1024), save :: ensemble_last_error = ''   ! facade-shared last-error text (get_last_bmi_error)

   real(real64), allocatable, save, target :: gwl(:)          ! MODFLOW head, metres (driver writes)
   real(real64), allocatable, save, target :: qbot_volume(:)  ! recharge depth over step, metres (we write)
   real(real64), allocatable, save, target :: storage_coef(:) ! specific yield, - (we write)

contains

   integer function ensemble_init(config_file, ncol_in) result(rc)
      character(len=*), intent(in) :: config_file
      integer,          intent(in) :: ncol_in
      type(error_collection_t) :: errors
      type(diag_overrides_t)   :: toml_ov
      integer :: i
      rc = 0
      call read_logging_overrides_from_file(config_file, toml_ov)
      call init_logging(default_embedded_config(), toml_ov)
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
         columns(i)%diag%instance_id = i
         columns(i)%soilwater%flcoupled_gwl = .true.
      end do

      allocate(gwl(ncol),          source=0.0_real64)
      allocate(qbot_volume(ncol),  source=0.0_real64)
      allocate(storage_coef(ncol), source=real(STORAGE_COEF_DEFAULT_X100, real64)/100.0_real64)
      first_step_done = .false.   ! re-seat columns to per-cell gwl on first step
      coupled_mode    = .true.

      if (.not. ensemble_all_swbotb1()) rc = 3
   end function ensemble_init

   !> Advance every column by exactly ONE day and report the per-day recharge
   !! depth (metres) into qbot_volume. This matches the MODFLOW coupler, which
   !! calls solve() once per DAILY stress period and reads qbot_volume as the
   !! recharge over that day (then divides by delt to get a rate).
   !!
   !! swap_run_step advances ONE Richards SUBSTEP (dt ~ 0.04 d here), not a full
   !! day: timecontrol_advance at the end of each step sets time%flDayEnd .true.
   !! only on the last substep of the day. So we loop swap_run_step until the
   !! day boundary and integrate the bottom flux over the substeps ourselves as
   !! sum(qbot*dt). dt MUST be read BEFORE swap_run_step: timecontrol_advance
   !! (at the end of the step) overwrites time%dt with the NEXT substep's dt,
   !! and qbot is the flux that applied over the dt that was current at entry.
   !! Integrating per-substep with each substep's own dt is exactly how
   !! waterbalance accumulates soil%cqbot; we replicate it locally so the result
   !! is immune to cqbot's flZeroCumu reset gating (cqbot is zeroed at the start
   !! of balance-output days, which would corrupt a naive cqbot delta).
   !!
   !! Sign convention: soil%qbot (cm/d) is NEGATIVE for downward percolation
   !! (water leaving the column into groundwater). We sum qbot*dt (negative for
   !! a recharging day) then flip the sign, so positive qbot_volume = recharge
   !! INTO groundwater, in metres.
   integer function ensemble_step_day() result(rc)
      integer      :: i
      real(real64) :: qbot_cm_day, dt_days
      rc = 0

      ! First coupled step: the driver has just injected each column's own
      ! per-cell water table into gwl(:) (metres). Every column was built at the
      ! shared config's uniform gwli, so re-seat its moisture profile to
      ! hydrostatic equilibrium with its actual gwl before solving. This removes
      ! the day-1 re-equilibration flux that otherwise seeds a coupled-solver
      ! oscillation (see reequilibrate_column_to_gwl). m -> cm.
      if (.not. first_step_done) then
         do i = 1, ncol
            call reequilibrate_column_to_gwl(columns(i), gwl(i) * 100.0_real64)
         end do
         first_step_done = .true.
      end if

      do i = 1, ncol
         ! Drive the whole day with the injected head (gwl_injected persists
         ! across substeps; BoundBottom reads it every substep).
         columns(i)%soilwater%gwl_injected = gwl(i) * 100.0_real64   ! m -> cm
         qbot_cm_day = 0.0_real64                                     ! cm, this day

         ! Advance exactly one day: run substeps until the day boundary
         ! (flDayEnd) or the simulation end (flRunEnd), integrating qbot*dt
         ! with the dt that was current for each substep (read before the step).
         do
            dt_days = columns(i)%timecontrol%dt
            call swap_run_step(columns(i), configs(column_config(i)))
            if (library_fatal_raised()) then; rc = 1; return; end if
            if (columns(i)%diag%aborted()) then; rc = 1; exit; end if
            qbot_cm_day = qbot_cm_day + columns(i)%soilwater%qbot * dt_days
            if (columns(i)%timecontrol%flDayEnd) exit
            if (columns(i)%timecontrol%flRunEnd) exit
         end do

         ! cm -> m, sign flip so positive qbot_volume = recharge into groundwater.
         qbot_volume(i) = -qbot_cm_day * 0.01_real64
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
      coupled_mode = .false.
   end function ensemble_finalize

   !> Allocate a 1-column, uncoupled ensemble (no exchange arrays, no gwl
   !! injection). The caller then populates configs(1) (from file or from an
   !! in-memory TOML string via the CAPI) and inits columns(1). This is the
   !! backing store for the single-column BMI/CAPI facade — the former
   !! capi_state singleton (ADR 0050 sub-arc 2).
   subroutine ensemble_allocate_single()
      ! Drop any stale storage without swap_close side effects (re-init in one
      ! process replaces the previous instance, matching the old singleton).
      if (allocated(columns))       deallocate(columns)
      if (allocated(configs))       deallocate(configs)
      if (allocated(column_config)) deallocate(column_config)
      allocate(configs(1))
      allocate(columns(1))
      allocate(column_config(1)); column_config = 1
      ncol = 1
      coupled_mode = .false.
   end subroutine ensemble_allocate_single

   !> Single-column standalone init from a TOML file — the merged facade's
   !! `initialize()` path when no ensemble.txt sidecar is present. Mirrors the
   !! former singleton path exactly: swap_init on one (state, config) pair.
   integer function ensemble_init_single(config_file) result(rc)
      character(len=*), intent(in) :: config_file
      rc = 0
      call ensemble_allocate_single()
      call swap_init(config_file, columns(1), configs(1))
   end function ensemble_init_single

   logical function ensemble_is_coupled()
      ensemble_is_coupled = coupled_mode
   end function ensemble_is_coupled

   function ensemble_column1() result(p)
      type(swap_state_t), pointer :: p
      p => null()
      if (allocated(columns)) p => columns(1)
   end function ensemble_column1

   function ensemble_config1() result(p)
      type(swap_config_t), pointer :: p
      p => null()
      if (allocated(configs)) p => configs(1)
   end function ensemble_config1

   !> Read the integer column count from the `ensemble.txt` sidecar next to the
   !! config (accepts a dir or a file path). Returns 0 when absent/invalid —
   !! the merged facade treats 0 as "single-column standalone mode".
   integer function read_ncol_sidecar(dir_or_file) result(n)
      character(len=*), intent(in) :: dir_or_file
      integer :: u, ios
      character(len=512) :: path
      n = 0
      path = trim(dir_or_file)
      if (index(path, '.toml') > 0) path = path(1:scan(path, '/', back=.true.))
      open (newunit=u, file=trim(path)//'ensemble.txt', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read (u, *, iostat=ios) n
      close (u)
      if (ios /= 0) n = 0
   end function read_ncol_sidecar

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
