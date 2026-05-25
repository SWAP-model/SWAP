! cropfixed_runtime.f90 — type-1 (fixed/simple) crop runtime dispatcher.
! ----------------------------------------------------------------------
      module cropfixed_runtime_mod
      implicit none
      private

      public :: cropfixed

      contains

! ----------------------------------------------------------------------
      subroutine cropfixed (task, state)
! ----------------------------------------------------------------------
!     date               : august 2004
!     purpose            : simple crop growth routine for swap
!
! [GR-CROP 2026-05-25] crop-sweep:
!   - magrs sourced from swap_array_dimensions.
!   - rdmax read via state%cfg%crop%rdmax (Class B direct read).
!   - siccaplai*lai branches: swinter=3 stub-errored on TOML path, siccaplai
!     always 0; replaced with 0.0d0 inline.
!   - mrftb/wrtb tables always zero on TOML path (legacy reader removed);
!     state%crop%oxygen%max_resp_factor / state%crop%oxygen%w_root_ss
!     directly assigned 0.
!   - Cross-file legacy globals still in use-list: wiltpoint/twilt/flhydrlift/reltr
!     are written here and read by rootextraction/cropgrass/cropwofost/jongvanlier.
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: magrs
      use variables, only: wiltpoint, twilt, flhydrlift,    &  ! cross-file with rootextraction/cropgrass/cropwofost/jongvanlier
                           reltr                              ! cross-file with cropgrass/cropwofost
      use soilhydraulics_utils, only: watcon
      use array_utils, only: afgen
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer   i,task,lcc,swhydrlift
      real(8)   dummy,dtsum,dvr

! --- rooting
      real(8)   rrpot,rr

      save
! ----------------------------------------------------------------------
      ! [GR-CROP 2026-05-25] sub-record associate (crop-sweep convention).
      associate( crop     => state%crop,            &
                 soil     => state%soilwater,       &
                 mesh     => state%mesh,            &
                 atmo     => state%atmosphere,      &
                 time     => state%timecontrol,     &
                 cfg_crop => state%cfg%crop         )

      select case (task)
      case (1)

! === initialization ===================================================
      
! --- read crop data: dispatch on per-rotation typed-config cache
!     (ADR 0016). Falls back to legacy reader for rotations whose
!     .crp.toml is not yet authored or for rotation types not yet
!     ported (Phase 1: only type=1 cropfixed; Phases 2/3 add types
!     2 and 3). Teardown: end of Phase 4 removes the else-branch.
      block
         use crop_config_global_mod, only: crop_config_global
         use cropfixed_init_mod, only: cropfixed_init_from_config
         logical :: use_cache
         use_cache = .false.
         if (associated(crop_config_global)) then
            if (allocated(crop_config_global%rotation_loaded)) then
               if (crop%common%icrop >= 1 .and. crop%common%icrop <= size(crop_config_global%rotation_loaded)) then
                  if (crop_config_global%rotation_loaded(crop%common%icrop)) use_cache = .true.
               end if
            end if
         end if
         if (use_cache) then
            call cropfixed_init_from_config(crop_config_global%rotation_fixed(crop%common%icrop), crop%common%icrop, lcc, state)
            ! swhydrlift is read by legacy readcropfixed only inside the
            ! swdrought=2 branch (stub-errored in Phase 1). Set to 0 here
            ! to mirror the default; Phase 2 (cropwofost) reuses this
            ! same field on its own dispatch path.
            swhydrlift = 0
         else
            ! ADR 0016 cache-miss: typed config required for type=1 rotations.
            ! No silent legacy fallback — the user must author cropfixed.crp.toml.
            call fatalerr_collected('cropgrowth/CropFixed', &
               'cropfixed rotation has no loaded .crp.toml — author the file or use the legacy executable.')
         end if
      end block

! --- maximum rooting depth
      if (crop%common%swrd.eq.1) then
        crop%common%rdm = cfg_crop%rdmax
      else
        crop%common%rdm = min(cfg_crop%rdmax,crop%common%rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (time%t1900 - time%tstart .gt. tiny .or. soil%swinco .ne. 3 .or.           &
     &  dabs(time%t1900 - crop%common%cropstart) .lt. tiny) then

        crop%common%dvs = 0.0d0

! --- actual rooting depth
        if (crop%common%swrd.eq.1) then
          crop%common%rd = afgen (crop%common%rdtb,22,crop%common%dvs)
          crop%common%rd = min(crop%common%rd,crop%common%rdm)
        else
          crop%common%rd = min(crop%common%rdi,crop%common%rdm)
        endif
        crop%common%rdpot = crop%common%rd

      endif

! --- initial lai or sc
      crop%lai = afgen (crop%fixed%gctb,(2*magrs),crop%common%dvs)
      if (crop%common%swgc.eq.2) then
        crop%common%gc = crop%lai
        crop%lai = crop%lai*3.0d0
      endif

! --- initial crop factor or crop height
      crop%common%cf = afgen (crop%fixed%cftb,(2*magrs),crop%common%dvs)
      crop%common%ch = afgen (crop%fixed%chtb,(2*magrs),crop%common%dvs)
      if (crop%swcf.eq.3) then
        crop%fixed%cfeic = afgen (crop%fixed%cfeictb,(2*magrs),crop%common%dvs)
      endif

! --- initial storage on canopy
      if (crop%common%swinter.eq.3) then
        atmo%siccapact = 0.0d0   ! [GR-CROP 2026-05-25] swinter=3 (Gash) stub-errored on TOML; siccaplai always 0
      endif

! --- initial dry weight of roots at soil surface; oxygen module
!     [GR-CROP 2026-05-25] wrtb always-zero on TOML (legacy reader removed) → 0
      state%crop%oxygen%w_root_ss = 0.0d0

! --- initial ratio root total respiration / maintenance respiration; oxygen module
!     [GR-CROP 2026-05-25] mrftb always-zero on TOML (legacy reader removed) → 0
      state%crop%oxygen%max_resp_factor = 0.0d0

! --- initialize matric flux potential (SS-CRP C-2.5: hroot/hleaf/mfluxtable
!     init moved to CropGrowth dispatcher which has access to state).
      if (crop%common%swdrought .eq. 2) then
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,mesh%numnod
         twilt(i) = watcon(wiltpoint, &
                            soil%vg_params(i), &
                            soil%iHWCKmodel(soil%layer(i)), &
                            i, soil)
        enddo
      endif

      return

      case (2)
      continue

! === calculate potential rate and state variables ======================
      case (3)

! === calculate actual rate and state variables ======================

! --- increase in temperature sum
      dtsum = max (0.0d0,atmo%Tav-crop%common%tbase)

! --- development rate
      if (crop%common%idev.eq.1) then
        dvr = 2.0/lcc
      elseif (crop%common%idev.eq.2) then
        if (crop%common%dvs.lt.1.0d0) then
          dvr = dtsum/crop%common%tsumea
        else
          dvr = dtsum/crop%common%tsumam
        endif
      endif

! --- water stress
      if(dabs(atmo%ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(min(soil%tra/atmo%ptra,1.0d0),0.0d0)
      endif

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      crop%common%dvs = min(crop%common%dvs+dvr,2.d0)
      crop%common%tsum = crop%common%tsum + dtsum

! --- leaf area index or soil cover fraction
      crop%lai = afgen (crop%fixed%gctb,(2*magrs),crop%common%dvs)
      if (crop%common%swgc.eq.2) then
        crop%common%gc = crop%lai
        crop%lai = crop%lai*3.0d0
      endif

! --- crop factor or crop height
      crop%common%cf        = afgen (crop%fixed%cftb,(2*magrs),crop%common%dvs)
      crop%common%ch        = afgen (crop%fixed%chtb,(2*magrs),crop%common%dvs)
      if (crop%swcf.eq.3) then
        crop%fixed%cfeic = afgen (crop%fixed%cfeictb,(2*magrs),crop%common%dvs)
      endif

! --- update canopy storage capacity
      if (crop%common%swinter.eq.3) then
        atmo%siccapact = 0.0d0   ! [GR-CROP 2026-05-25] swinter=3 (Gash) stub-errored on TOML; siccaplai always 0
      endif

! --- dry weight of roots at soil surface; oxygen module
!     [GR-CROP 2026-05-25] wrtb always-zero on TOML (legacy reader removed) → 0
      state%crop%oxygen%w_root_ss = 0.0d0

! --- ratio root total respiration / maintenance respiration; oxygen module
!     [GR-CROP 2026-05-25] mrftb always-zero on TOML (legacy reader removed) → 0
      state%crop%oxygen%max_resp_factor = 0.0d0

      case (4)

! --- root extension
      if (crop%common%swrd.eq.1) then
        crop%common%rdpot = afgen (crop%common%rdtb,22,crop%common%dvs)
        crop%common%rdpot = min(crop%common%rdpot,crop%common%rdm)
        crop%common%rd    = crop%common%rdpot
      else
        rrpot = min (crop%common%rdm-crop%common%rdpot,crop%common%rri)
        if (atmo%ptra.lt.nihil) rrpot = 0.0d0
        crop%common%rdpot = crop%common%rdpot + rrpot

        rr = min (crop%common%rdm-crop%common%rd,crop%common%rri)
        if (atmo%ptra.lt.nihil .or. soil%flWrtNonox) rr = 0.0d0
        if (crop%common%swdmi2rd.eq.1 .and. atmo%ptra.ge.nihil) rr = rr * soil%tra/atmo%ptra
        crop%common%rd = crop%common%rd + rr
      endif

      return

      case default
         call fatalerr_collected ('CropFixed', 'Illegal value for TASK')
      end select

      end associate
      return
      end subroutine cropfixed

      end module cropfixed_runtime_mod
