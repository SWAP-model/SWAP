!> @file crop_state.f90
!! SS-GR-ATM: typed crop runtime state. Foundation for Arc 8 (crop
!! cluster). Arc 4 (atmosphere) introduces this module to host the
!! ~12 crop-runtime symbols read by atmosphere code; Arc 8 later
!! migrates all remaining crop readers.
module crop_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use crop_common_state_mod,      only: crop_common_state_t
   use crop_fixed_state_mod,       only: crop_fixed_state_t
   use crop_wofost_state_mod,      only: crop_wofost_state_t
   use crop_grass_state_mod,       only: crop_grass_state_t
   use crop_irrigation_state_mod,  only: crop_irrigation_state_t
   use crop_oxygen_state_mod,      only: crop_oxygen_state_t
   use cropfixed_config_mod,       only: cropfixed_config_t
   use cropgrass_config_mod,       only: cropgrass_config_t
   use cropwofost_config_mod,      only: cropwofost_config_t
   implicit none
   private
   public :: crop_state_t

   type :: crop_state_t
      real(real64) :: lai             = 0.0_real64    !! leaf area index (-)
      real(real64) :: kdif            = 0.0_real64    !! extinction for diffuse light (-)
      real(real64) :: kdir            = 0.0_real64    !! extinction for direct light (-)
      real(real64) :: cofab           = 0.0_real64    !! interception coefficient (cm)
      real(real64) :: cfbs            = 1.0_real64    !! bare-soil ET factor (-)
      integer      :: swcf            = 0             !! crop factor switch
      integer      :: swcfbs          = 0             !! bare-soil factor switch
      real(real64) :: gird            = 0.0_real64    !! gross irrigation depth (cm)
      logical      :: flCropEmergence = .false.       !! crop emerged flag
      real(real64) :: et0             = 0.0_real64    !! potential ET (cm/d)
      real(real64) :: ew0             = 0.0_real64    !! potential evap wet crop (cm/d)
      real(real64) :: es0             = 0.0_real64    !! potential evap bare soil (cm/d)
      ! [SS-GR-CROP A2] sub-record for shared crop runtime fields
      type(crop_common_state_t) :: common
      ! [SS-GR-CROP A3] sub-record for fixed-crop runtime fields
      type(crop_fixed_state_t) :: fixed
      ! [SS-GR-CROP A4] sub-record for WOFOST biomass pools + flows
      type(crop_wofost_state_t) :: wofost
      ! [SS-GR-CROP A5] sub-record for grass-specific runtime fields
      type(crop_grass_state_t) :: grass
      ! [GR-CROP 2026-05-25] SSDI persistent state — migrated from variables.f90 _irr globals
      type(crop_irrigation_state_t) :: irrigation
      ! [GR-CROP 2026-05-25] Bartholomeus oxygen-stress workspace + SAVE state
      ! migrated from O2_pars module + o2_* legacy globals
      type(crop_oxygen_state_t) :: oxygen
      ! [state%cfg-retirement cluster 6] per-rotation metadata; snapshotted in crop_state_init.
      real(real64), allocatable :: rotation_start(:)  !! crop season start (days since 1900), per rotation; snapshotted from config%crop%rotation_start
      real(real64), allocatable :: rotation_end(:)    !! crop season end (days since 1900), per rotation; snapshotted from config%crop%rotation_end
      logical,      allocatable :: rotation_loaded(:) !! per-rotation typed-config cache available flag; snapshotted from config%crop%rotation_loaded

      ! [state%cfg-retirement mop-up] per-rotation typed-config sub-arrays; snapshotted in crop_state_init.
      ! Compute paths in cropgrowth_helpers / cropwofost_runtime / cropgrass_runtime /
      ! cropfixed_runtime / cropgrowth read these via state%crop instead of state%cfg%crop.
      type(cropfixed_config_t),  allocatable :: rotation_fixed(:)   !! cropfixed parameters per rotation
      type(cropgrass_config_t),  allocatable :: rotation_grass(:)   !! cropgrass parameters per rotation
      type(cropwofost_config_t), allocatable :: rotation_wofost(:)  !! cropwofost parameters per rotation
   contains
      procedure :: init => crop_state_init
   end type crop_state_t

contains

   subroutine crop_state_init(self, crop_cfg, meteo_cfg, pathwork_in)
      use crop_config_mod,              only: crop_config_t
      use meteorology_config_mod,       only: meteorology_config_t
      class(crop_state_t),         intent(inout) :: self
      type(crop_config_t), target, intent(in)    :: crop_cfg
      type(meteorology_config_t),  intent(in)    :: meteo_cfg
      character(len=*),            intent(in)    :: pathwork_in
      integer :: i

      ! [GR-SEED 2026-05-25 Task 10] Absorbed from config_to_variables adapter:

      ! 1. swcrop==1 flags: arm per-rotation reader gates so cropgrowth.f90's
      !    per-crop init (ArableLandGerm/CropFixed/Wofost/Grass) actually fires.
      !    Without this, flCropReadFile stays .false. and every rotation is treated
      !    as bare soil: LAI/cf/rd remain 0, TPOT/TACT collapse.
      if (crop_cfg%swcrop == 1) then
         self%common%flCropReadFile = .true.
         self%common%flCropOpenFile = .true.
      end if

      ! 2. rotation_type → state%crop%common%croptype
      if (allocated(crop_cfg%rotation_type)) then
         if (.not. allocated(self%common%croptype)) then
            allocate(self%common%croptype(size(crop_cfg%rotation_type)))
         end if
         do i = 1, size(crop_cfg%rotation_type)
            self%common%croptype(i) = crop_cfg%rotation_type(i)
         end do
      end if

      ! 3. Evaporation cfbs (deferred from Task 2)
      self%cfbs = meteo_cfg%evaporation%cfbs

      ! 4. rdmax — single config scalar consumed by cropgrass/cropfixed init (task=1)
      self%common%rdmax = crop_cfg%rdmax

      ! 5. rotation_start / rotation_end — per-rotation time windows (days since 1900).
      !    Used by cropgrowth.f90 task=1 to find the active rotation and task=2
      !    for harvest detection.
      if (allocated(crop_cfg%rotation_start)) then
         if (.not. allocated(self%rotation_start)) then
            allocate(self%rotation_start(size(crop_cfg%rotation_start)))
         end if
         self%rotation_start(:) = crop_cfg%rotation_start(:)
      end if
      if (allocated(crop_cfg%rotation_end)) then
         if (.not. allocated(self%rotation_end)) then
            allocate(self%rotation_end(size(crop_cfg%rotation_end)))
         end if
         self%rotation_end(:) = crop_cfg%rotation_end(:)
      end if
      if (allocated(crop_cfg%rotation_loaded)) then
         if (.not. allocated(self%rotation_loaded)) then
            allocate(self%rotation_loaded(size(crop_cfg%rotation_loaded)))
         end if
         self%rotation_loaded(:) = crop_cfg%rotation_loaded(:)
      end if

      ! [state%cfg-retirement mop-up] Snapshot typed rotation sub-arrays so compute
      ! paths can read state%crop%rotation_{wofost,fixed,grass} without touching state%cfg.
      if (allocated(crop_cfg%rotation_wofost)) then
         if (.not. allocated(self%rotation_wofost)) then
            allocate(self%rotation_wofost(size(crop_cfg%rotation_wofost)))
         end if
         self%rotation_wofost(:) = crop_cfg%rotation_wofost(:)
      end if
      if (allocated(crop_cfg%rotation_fixed)) then
         if (.not. allocated(self%rotation_fixed)) then
            allocate(self%rotation_fixed(size(crop_cfg%rotation_fixed)))
         end if
         self%rotation_fixed(:) = crop_cfg%rotation_fixed(:)
      end if
      if (allocated(crop_cfg%rotation_grass)) then
         if (.not. allocated(self%rotation_grass)) then
            allocate(self%rotation_grass(size(crop_cfg%rotation_grass)))
         end if
         self%rotation_grass(:) = crop_cfg%rotation_grass(:)
      end if

      ! pathwork_in reserved for future CSV seed migration (consistency with sibling inits).
      ! Not consumed yet; gfortran does not warn on unused dummy args by default.
   end subroutine crop_state_init

end module crop_state_mod
