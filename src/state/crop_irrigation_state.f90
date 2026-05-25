!> @file crop_irrigation_state.f90
!! SS-GR-CROP: typed crop runtime state — irrigation runtime/persistent state.
!! Populated by state%crop%irrigation%init (called from swap_mod after CalcGrid);
!! mutated at runtime by src/crop/irrigation.f90 SSDI_irrigation(2)/SSDI_irrigation(9).
!!
!! [GR-CROP 2026-05-25] hosts the 16 SSDI persistent state fields that
!! were SAVE-state on irrigation.f90 module variables (later legacy
!! globals in variables.f90 module). The corresponding `*_irr` legacy
!! globals retire in the same commit.
!!
!! dt_SSDI_event lives on this sub-record as of 2026-05-25 — the
!! cross-file consumer src/core/timecontrol_mod.f90 was migrated as
!! part of the irrigation.f90 sub-arc follow-up.
!!
!! [GR-SEED 2026-05-25 Task 5] init signature extended to absorb:
!!   - fixed-irrigation event seeding (seed_fixed_irrigation, private)
!!   - SSDI seeding (apply_ssdi_seed, public test seam;
!!                   apply_ssdi_mode0/apply_ssdi_mode1, private)
module crop_irrigation_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: mairg
   implicit none
   private
   public :: crop_irrigation_state_t
   public :: apply_ssdi_seed

   type :: crop_irrigation_state_t

      ! SSDI configuration (init-once from state%cfg%irrigation%ssdi via
      ! state%crop%irrigation%init called after CalcGrid in swap_mod)
      integer      :: swssdi            = 0          !! SSDI active (0=no, 1=yes)
      integer      :: nod_ssdi(2)       = 0          !! Upper and lower nodes for SSDI
      integer      :: ssdi_schedule     = 0          !! Schedule type (0=fixed dates, 1=internal)
      integer      :: ssdi_sched_type   = 0          !! Internal schedule type (1=Tact/Tpot, 2=h, 3=theta)
      integer      :: nod_ssdi_sensor   = 0          !! Sensor node (if ssdi_sched_type > 1)
      real(real64) :: ssdi_threshold    = 0.0_real64 !! Threshold value for scheduling
      real(real64) :: ssdi_threshold_z  = 0.0_real64 !! Depth for threshold value (cm)
      real(real64) :: ssdi_amount       = 0.0_real64 !! Amount of scheduled irrigation (cm)
      real(real64) :: ssdi_appl_rate    = 0.0_real64 !! Application rate (cm/d)
      integer      :: sw_interval       = 0          !! Switch for minimum interval
      integer      :: days_interval     = 1          !! Minimum days between applications

      ! SSDI runtime cursor/counter state
      integer      :: days_counter      = 366        !! Days since previous application
      integer      :: nirri             = 1          !! SSDI counter / entry point into ssdi_date
      real(real64) :: dt_SSDI_event     = 1.0_real64 !! Sub-day SSDI event end fraction (1.0 = no event)

      ! Fixed-date schedule tables (mode 0)
      real(real64) :: ssdi_date(mairg)     = 0.0_real64 !! Fixed irrigation dates (days-since-1900)
      real(real64) :: ssdi_rate_f(mairg)   = 0.0_real64 !! Fixed irrigation rates (cm/d)
      real(real64) :: ssdi_amount_f(mairg) = 0.0_real64 !! Fixed irrigation amounts (cm)

      ! Surface fixed-irrigation events (populated by seed_fixed_irrigation
      ! from config%irrigation%fixed_events / fixed_events_file when swirfix == 1).
      integer      :: nirri_fixed              = 1          !! Cursor into fixed-irrigation event arrays
      real(real64) :: irdate(mairg)            = 0.0_real64 !! Fixed irrigation dates (days-since-1900)
      real(real64) :: irdepth(mairg)           = 0.0_real64 !! Fixed irrigation depths (cm)
      real(real64) :: irconc(mairg)            = 0.0_real64 !! Fixed irrigation concentrations (M/L3)
      integer      :: irtype(mairg)            = 0          !! Fixed irrigation types (0=sprinkler, 1=surface)

   contains
      procedure :: init => crop_irrigation_state_init
   end type crop_irrigation_state_t

contains

   !> Extended init: seed fixed-irrigation events + SSDI persistent state
   !! from config. Called after CalcGrid so mesh node-resolution is valid.
   subroutine crop_irrigation_state_init(self, config_irrigation, tstart, tend, mesh, pathwork_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_config_t
      use mesh_state_mod,        only: mesh_state_t
      class(crop_irrigation_state_t), intent(inout) :: self
      type(irrigation_config_t),      intent(in)    :: config_irrigation
      real(real64),                   intent(in)    :: tstart, tend
      type(mesh_state_t),             intent(in)    :: mesh
      character(len=*),               intent(in)    :: pathwork_in

      ! Defaults are already set on the type declaration; mirror top-level flags.
      self%swssdi = config_irrigation%swssdi
      self%nirri  = 1

      ! Fixed-irrigation events seeding (formerly adapter Irrigation block).
      call seed_fixed_irrigation(self, config_irrigation, pathwork_in)

      ! SSDI seeding (formerly apply_irrigation_ssdi in config_to_variables.f90).
      if (config_irrigation%swssdi == 1) then
         call apply_ssdi_seed(self, config_irrigation%ssdi, tstart, tend, mesh, pathwork_in)
      end if
   end subroutine crop_irrigation_state_init


   !> Seed fixed-irrigation event arrays from config inline table or CSV file.
   !! Formerly the "Irrigation (audit: 8 fields)" block in config_to_variables.f90.
   !! mm->cm conversion on irdepth mirrors readswap.f90:503.
   subroutine seed_fixed_irrigation(self, config_irrigation, pathwork_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_config_t
      type(crop_irrigation_state_t), intent(inout) :: self
      type(irrigation_config_t),     intent(in)    :: config_irrigation
      character(len=*),              intent(in)    :: pathwork_in
      integer :: i, n

      ! Inline fixed events: copy (date, depth, conc, type) rows from
      ! the typed config table into the state arrays. The reader already
      ! stored col 1 as days-since-1900, so this is a straight copy.
      ! The /10.0 on irdepth mirrors readswap.f90:503 (mm -> cm).
      if (allocated(config_irrigation%fixed_events)) then
         n = size(config_irrigation%fixed_events, 1)
         do i = 1, min(n, size(self%irdate))
            self%irdate(i)  = config_irrigation%fixed_events(i, 1)
            self%irdepth(i) = config_irrigation%fixed_events(i, 2) / 10.0d0
            self%irconc(i)  = config_irrigation%fixed_events(i, 3)
            self%irtype(i)  = nint(config_irrigation%fixed_events(i, 4))
         end do
      else if (config_irrigation%swirfix == 1 .and. allocated(config_irrigation%fixed_events_file)) then
         ! Phase 4f cleanup: long-form fixed-irrigation events outsourced
         ! to a CSV companion file (date, depth_mm, conc, type). The
         ! reader emits days-since-1900 in column 1; the rest of the
         ! unpack mirrors the inline-fixed_events path above (mm -> cm
         ! on depth, nint() on type). Replaces the legacy .irg HACK.
         if (len_trim(config_irrigation%fixed_events_file) > 0) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               use iso_fortran_env, only: real64
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t) :: csv_errs
               integer :: k_csv, nrows_csv
               character(len=5) :: irrig_header(4)
               irrig_header(1) = 'date '
               irrig_header(2) = 'depth'
               irrig_header(3) = 'conc '
               irrig_header(4) = 'type '
               call read_csv_table( &
                  trim(pathwork_in)//trim(config_irrigation%fixed_events_file), &
                  irrig_header, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows_csv = size(csv_table, 1)
                  do k_csv = 1, min(nrows_csv, size(self%irdate))
                     self%irdate(k_csv)  = csv_table(k_csv, 1)
                     self%irdepth(k_csv) = csv_table(k_csv, 2) / 10.0d0  ! mm -> cm
                     self%irconc(k_csv)  = csv_table(k_csv, 3)
                     self%irtype(k_csv)  = nint(csv_table(k_csv, 4))
                  end do
               end if
            end block
         end if
      end if
      ! [GR-CROP 2026-05-25] nirri_fixed canonical cursor. Type default = 1;
      ! reaffirmed here so swap_mod re-entries get a consistent reset.
      self%nirri_fixed = 1
   end subroutine seed_fixed_irrigation


   !> Seed SSDI persistent state from config. Resolves ssdi_z depths to mesh
   !! node indices, then dispatches to mode0 (fixed-date CSV) or mode1
   !! (scheduled-trigger). Public: used as a test seam by
   !! tests/unit/io/toml/test_apply_irrigation_ssdi.pf.
   !!
   !! Formerly apply_irrigation_ssdi in config_to_variables.f90.
   subroutine apply_ssdi_seed(self, ssdi, tstart, tend, mesh, pathwork_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_ssdi_t
      use mesh_state_mod,        only: mesh_state_t
      use error_mod,             only: fatalerr_collected
      type(crop_irrigation_state_t), intent(inout) :: self
      type(irrigation_ssdi_t),       intent(in)    :: ssdi
      real(real64),                  intent(in)    :: tstart, tend
      type(mesh_state_t),            intent(in)    :: mesh
      character(len=*),              intent(in), optional :: pathwork_in

      integer :: i, j, nod_top, nod_bot, ncomp

      ! Resolve ssdi_z(1:2) -> layer indices via zbotcp walk.
      ! Mirrors the legacy SSDI_irrigation(1) loop at irrigation.f90:387-393.
      ! [GR-BH Task 35] zbotcp/NumNod globals replaced by mesh fields.
      ! Guard: mesh not yet populated at config_to_variables call time;
      ! self%nod_ssdi defaults to 0 if mesh not built (resolved after CalcGrid).
      self%nod_ssdi = 0
      if (mesh%numnod > 0 .and. allocated(mesh%zbotcp)) then
         do j = 1, 2
            i = 1
            do while (mesh%zbotcp(i) > (ssdi%ssdi_z(j) + 1.0e-5_real64))
               i = i + 1
               if (i > mesh%numnod) exit
            end do
            self%nod_ssdi(j) = i
         end do
      end if
      nod_top = self%nod_ssdi(1)
      nod_bot = self%nod_ssdi(2)
      ncomp   = nod_bot - nod_top + 1
      if (ncomp < 1) then
         call fatalerr_collected('apply_ssdi_seed', &
                                 'ssdi_z resolves to zero compartments — check ssdi_z vs grid')
      end if

      ! [GR-SOIL 2026-05-24] qssdi zero-init handled by soilwater_init (state field).
      self%dt_SSDI_event = 1.0_real64

      select case (ssdi%schedule)
      case (0)
         if (present(pathwork_in)) then
            call apply_ssdi_mode0(ssdi, ncomp, tstart, tend, self, pathwork_in)
         else
            call apply_ssdi_mode0(ssdi, ncomp, tstart, tend, self, './')
         end if
      case (1)
         call apply_ssdi_mode1(ssdi, ncomp, self)
      end select
   end subroutine apply_ssdi_seed


   !> Mode-0 (fixed-date): stage CSV; populate ssdi_*_f_irr; deferred
   !! date-window validation; initial nirri_ssdi_irr entry-point from tstart.
   !! Formerly apply_ssdi_mode0 in config_to_variables.f90.
   subroutine apply_ssdi_mode0(ssdi, ncomp, tstart, tend, self, pathwork_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use csv_reader_mod, only: read_csv_table
      use error_mod, only: error_collection_t, fatalerr_collected
      use irrigation_config_mod, only: irrigation_ssdi_t
      ! [GR-IO 2026-05-25 Phase 4] mairg from canonical swap_array_dimensions module
      use swap_array_dimensions, only: mairg
      type(irrigation_ssdi_t),       intent(in)    :: ssdi
      integer,                       intent(in)    :: ncomp
      real(real64),                  intent(in)    :: tstart, tend
      type(crop_irrigation_state_t), intent(inout) :: self
      character(len=*),              intent(in)    :: pathwork_in

      real(real64), allocatable :: tbl(:,:)
      type(error_collection_t)  :: errs
      character(len=300) :: csvpath
      character(len=8)   :: hdr(3)
      integer :: i, n, nirri_init
      logical :: any_in_window, window_in_dates

      hdr(1) = 'date    '
      hdr(2) = 'rate_f  '
      hdr(3) = 'amount_f'
      csvpath = trim(pathwork_in) // trim(ssdi%fixed%events_file)
      call read_csv_table(trim(csvpath), hdr, tbl, errs)
      call errs%abort_if_fatal()

      n = 0
      if (allocated(tbl)) n = size(tbl, 1)
      if (n < 1) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: events CSV has no rows')
         return
      end if
      if (n > mairg) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: events CSV exceeds mairg rows')
         return
      end if

      ! self%ssdi_date / ssdi_rate_f / ssdi_amount_f are
      ! fixed-size arrays of size(mairg); zero them and fill from the CSV.
      self%ssdi_date     = 0.0_real64
      self%ssdi_rate_f   = 0.0_real64
      self%ssdi_amount_f = 0.0_real64

      do i = 1, n
         self%ssdi_date(i) = tbl(i, 1)
         if (i > 1 .and. self%ssdi_date(i) <= self%ssdi_date(i-1)) then
            call fatalerr_collected('apply_irrigation_ssdi', &
                                    'mode 0: ssdi_date not strictly ascending')
            return
         end if
         ! mm/h -> cm/d (mirrors irrigation.f90:514)
         self%ssdi_rate_f(i)   = tbl(i, 2) * 0.1_real64 * 24.0_real64
         ! mm -> cm, then spread over `ncomp` compartments (mirrors irrigation.f90:397)
         self%ssdi_amount_f(i) = (tbl(i, 3) * 0.1_real64) / real(ncomp, real64)
      end do

      ! Date-window check: at least one date in [tstart, tend], OR
      ! [tstart, tend] contained in [date(1), date(n)].
      ! Replaces the deleted checkdate call (irrigation.f90:552).
      any_in_window  = .false.
      do i = 1, n
         if (self%ssdi_date(i) >= tstart - 1.0e-6_real64 .and. &
             self%ssdi_date(i) <= tend   + 1.0e-6_real64) then
            any_in_window = .true.
            exit
         end if
      end do
      window_in_dates = (self%ssdi_date(1) <= tstart + 1.0e-6_real64 .and. &
                         self%ssdi_date(n) >= tend   - 1.0e-6_real64)
      if (.not. any_in_window .and. .not. window_in_dates) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: no ssdi_date within simulation period')
      end if

      ! Determine initial entry point (mirrors irrigation.f90:368-373).
      nirri_init = 1
      do i = 1, n - 1
         if (tstart >= self%ssdi_date(i)) nirri_init = i
      end do
      if (tstart >= self%ssdi_date(n)) nirri_init = n
      self%nirri = nirri_init
   end subroutine apply_ssdi_mode0


   !> Mode-1 (scheduled-trigger): copy scheduled sub-block to state fields;
   !! set nirri/dt_SSDI_event/days_counter defaults
   !! (preserves the d8a88d6 regression-fix invariant for dt_SSDI_event = 1.0).
   !! Formerly apply_ssdi_mode1 in config_to_variables.f90.
   subroutine apply_ssdi_mode1(ssdi, ncomp, self)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_ssdi_t
      type(irrigation_ssdi_t),       intent(in)    :: ssdi
      integer,                       intent(in)    :: ncomp
      type(crop_irrigation_state_t), intent(inout) :: self

      self%ssdi_sched_type  = ssdi%scheduled%sched_type
      self%ssdi_threshold   = ssdi%scheduled%threshold
      self%ssdi_threshold_z = ssdi%scheduled%threshold_depth

      ! mm -> cm, then spread over ncomp compartments (mirrors irrigation.f90:399)
      self%ssdi_amount      = (ssdi%scheduled%ssdi_amount * 0.1_real64) / &
                               real(ncomp, real64)
      ! mm/h -> cm/d (mirrors irrigation.f90:546)
      self%ssdi_appl_rate   = ssdi%scheduled%ssdi_appl_rate * 0.1_real64 * 24.0_real64

      self%sw_interval      = ssdi%scheduled%sw_interval
      ! mirrors irrigation.f90:535-538
      if (ssdi%scheduled%sw_interval == 0) then
         self%days_interval = 1
      else
         self%days_interval = ssdi%scheduled%days_interval
      end if
      self%days_counter = 366   ! mirrors irrigation.f90:540

      self%nirri = 1
   end subroutine apply_ssdi_mode1

end module crop_irrigation_state_mod
