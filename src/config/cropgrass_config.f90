!> Type 3 (grass / WOFOST grass) crop config — populated from a .crp.toml file.
!! Field set is similar to cropfixed_config_t plus grass-specific management
!! (mowing, grazing, fertilizer). Per design Q2A, no shared base type.
module cropgrass_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   use irrigation_config_mod, only: irrigation_schedule_t
   implicit none
   private

   public :: cropgrass_config_t

   type :: cropgrass_config_t
      ! Phenology / development
      integer      :: idev = 2                  !! For grass typically 2 (temp-sum-based)
      integer      :: lcc  = 0
      real(real64) :: tbase = 0.0_real64
      real(real64) :: tsum1 = 0.0_real64
      real(real64) :: tsum2 = 0.0_real64

      ! Light & growth
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64) :: eff  = 0.0_real64
      real(real64) :: amax = 0.0_real64

      ! Tables
      real(real64), allocatable :: cftb(:)
      real(real64), allocatable :: chtb(:)
      real(real64), allocatable :: rdctb(:)

      ! Root growth
      real(real64) :: rdi = 0.0_real64
      real(real64) :: rri = 0.0_real64
      real(real64) :: rdc = 0.0_real64

      ! Water stress (Feddes) — same as fixed
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64

      ! Salinity & interception
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64
      real(real64) :: cofab  = 0.0_real64

      ! Mowing schedule (grass-specific)
      ! swharv: 0 = no scheduled mow, 1 = DM-threshold-driven, 2 = fixed-date table
      integer :: swharv = 0
      integer :: nmow   = 0                     !! Number of mowing events
      real(real64), allocatable :: dates_mowing(:)
      real(real64), allocatable :: lai_after_mow(:)

      ! Phase 4d additions for per-event mowing/grazing tables.
      ! Mowing block (when swharv=1 or 2):
      integer      :: swdmmow = 0                           !! 0=use heights, 1=DM threshold, 2=DM threshold (legacy SWDMMOW=2 in cases)
      real(real64), allocatable :: mowing_dates(:)          !! day-of-year per event (when swharv=2)
      real(real64), allocatable :: mowing_heights(:)        !! optional (when swdmmow=0)
      real(real64) :: dmharvest      = 0.0_real64           !! DM threshold (when swdmmow=1)
      real(real64) :: daylastharvest = 0.0_real64
      real(real64) :: dmlastharvest  = 0.0_real64
      integer      :: maxdaymow = 0

      ! Grazing (grass-specific)
      integer :: swgraz = 0                     !! 0=no grazing, 1=scheduled
      integer :: nstart_graz = 0                !! Start day-of-year
      integer :: nstop_graz  = 0                !! Stop day-of-year

      ! Phase 4d additions for grazing per-event details (when swgraz=1):
      integer      :: maxdaygrz = 0
      real(real64) :: dmgrazing = 0.0_real64
      integer      :: swdmgrz   = 0
      real(real64), allocatable :: lsdb(:)                  !! per-day stocking density
      real(real64) :: tagprest = 0.0_real64                 !! threshold above-ground residue

      ! Per-crop irrigation schedule (Phase 4d Task 12)
      type(irrigation_schedule_t) :: schedule
   contains
      procedure :: validate => cropgrass_config_validate
      procedure :: finalize => cropgrass_config_finalize
   end type cropgrass_config_t

contains

   subroutine cropgrass_config_validate(self, errors)
      class(cropgrass_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      call check_int_enum(self%idev,   [1, 2],       'cropgrass.idev',   errors)
      call check_int_enum(self%swharv, [0, 1, 2],    'cropgrass.swharv', errors)
      call check_int_enum(self%swgraz, [0, 1],       'cropgrass.swgraz', errors)

      call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'cropgrass.rdi', errors)
      call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'cropgrass.rdc', errors)
      call check_nonnegative_real(self%rri, 'cropgrass.rri', errors)

      call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropgrass', errors)
      call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropgrass', errors)

      call check_nonnegative_real(self%ecmax,  'cropgrass.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropgrass.ecslop', errors)

      ! ----- Mowing branches (swharv = 0/1/2) -----
      if (self%swharv == 1) then
         ! DM-threshold mowing: nmow drives event budget; mowing_dates not required
         call check_int_range(self%nmow, 1, 50, 'cropgrass.nmow', errors)
         call check_int_enum(self%swdmmow, [0, 1, 2], 'cropgrass.swdmmow', errors)
         if (self%swdmmow == 1) then
            if (self%dmharvest <= 0.0_real64) then
               call check_real_range(self%dmharvest, tiny(1.0_real64), huge(1.0_real64), &
                                     'cropgrass.dmharvest', errors)
            end if
         end if
         call check_nonnegative_real(self%daylastharvest, 'cropgrass.daylastharvest', errors)
         call check_nonnegative_real(self%dmlastharvest,  'cropgrass.dmlastharvest',  errors)
         call check_int_range(self%maxdaymow, 1, 366, 'cropgrass.maxdaymow', errors)
      else if (self%swharv == 2) then
         ! Fixed-date mowing: require mowing_dates allocated and nmow == size
         call check_int_range(self%nmow, 1, 50, 'cropgrass.nmow', errors)
         call check_int_enum(self%swdmmow, [0, 1, 2], 'cropgrass.swdmmow', errors)
         if (.not. allocated(self%mowing_dates)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                               'required when swharv=2 (allocate event table)', &
                               'cropgrass.mowing_dates')
         else if (size(self%mowing_dates) /= self%nmow) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                               'size must equal nmow', &
                               'cropgrass.mowing_dates')
         end if
         if (self%swdmmow == 0) then
            if (.not. allocated(self%mowing_heights)) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                  'required when swharv=2 and swdmmow=0', &
                                  'cropgrass.mowing_heights')
            else if (size(self%mowing_heights) /= self%nmow) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                  'size must equal nmow', &
                                  'cropgrass.mowing_heights')
            end if
         end if
      end if

      ! ----- Grazing branch (swgraz = 1) -----
      if (self%swgraz == 1) then
         call check_int_range(self%nstart_graz, 1, 366, 'cropgrass.nstart_graz', errors)
         call check_int_range(self%nstop_graz,  1, 366, 'cropgrass.nstop_graz',  errors)
         if (self%nstop_graz < self%nstart_graz) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                               'must be >= nstart_graz', &
                               'cropgrass.nstop_graz')
         end if
         call check_int_range(self%maxdaygrz, 1, 366, 'cropgrass.maxdaygrz', errors)
         if (self%dmgrazing <= 0.0_real64) then
            call check_real_range(self%dmgrazing, tiny(1.0_real64), huge(1.0_real64), &
                                  'cropgrass.dmgrazing', errors)
         end if
         call check_int_enum(self%swdmgrz, [0, 1, 2], 'cropgrass.swdmgrz', errors)
         call check_nonnegative_real(self%tagprest, 'cropgrass.tagprest', errors)
      end if

      ! Optional lsdb: if allocated, sanity-check non-negative entries
      if (allocated(self%lsdb)) then
         block
            integer :: i
            do i = 1, size(self%lsdb)
               if (self%lsdb(i) < 0.0_real64) then
                  call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                     'stocking density must be non-negative', &
                                     'cropgrass.lsdb')
                  exit
               end if
            end do
         end block
      end if

      call self%schedule%validate(errors)
   end subroutine cropgrass_config_validate

   subroutine cropgrass_config_finalize(self, errors)
      class(cropgrass_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      call self%schedule%finalize(errors)
   end subroutine cropgrass_config_finalize

end module cropgrass_config_mod
