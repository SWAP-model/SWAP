! cropfixed_runtime.f90
! GR-CROPWS Phase 0 Commit 0.4: cropfixed extracted from cropgrowth.f90.
! Pure relocation — no behavior change.
! ----------------------------------------------------------------------
      module cropfixed_runtime_mod
      implicit none
      private

      public :: cropfixed

      contains

!
! ----------------------------------------------------------------------
      subroutine cropfixed (task, state)
! ----------------------------------------------------------------------
!     date               : august 2004
!     purpose            : simple crop growth routine for swap
! SS-CRP C-2.5: state added (optional, intent in) to read flWrtNonox.
! SS-TC TC-10: t1900 read via state%timecontrol tc_t1900 alias.
! SS-GR-ATM A5.1: intent changed inout to allow dual-write in cropfixed_init_from_config.
! [SS-GR-CROPWS A2]: state optional removed — all callers pass state; all if(present(state)) guards dropped.
! [GR-CROP Phase B/6] narrow use variables
! [GR-CROPWS B1]: reads migrated to state%crop%X — icrop, dvs, tsum, rd, rdpot, rdm,
!   rdi, rri, rdc, ch, cf, cropstart, lai, swcf, cftb, chtb, cfeic, cfeictb.
!   Remaining in variables (no state home): magrs, idev, max_resp_factor, swrd, swgc,
!   swcf (keep for write), swinter, swdrought, swdmi2rd, tbase, tsumea, tsumam, rdmax,
!   siccapact, siccaplai, W_root_ss, wiltpoint, twilt, flhydrlift, gc, cfeic (write),
!   gctb, rdtb, mrftb, wrtb, swinco, reltr.
! ----------------------------------------------------------------------
      use variables, only: magrs, idev, &         ! [GR-CROPWS B1] reads→state; writes remain legacy; dvs/tsum/lai retired
                           max_resp_factor,                        &  ! rd/rdpot/swrd retired
                           swgc, swcf, swinter, swdrought,         &  ! swdmi2rd retired
                           tbase, tsumea, tsumam, rdmax,                  &
                           siccaplai, w_root_ss, wiltpoint,   &
                           twilt, flhydrlift, gc, cfeic,                 &
                           gctb, rdtb, mrftb, wrtb,                      &
                           swinco, reltr
      use soilhydraulics_utils, only: watcon
      use array_utils, only: afgen
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state   ! [SS-GR-CROPWS A2] removed optional — all callers pass state

! --- local variables
      integer   i,task,lcc,swhydrlift
      real(8)   dummy,dtsum,dvr

! --- rooting
      real(8)   rrpot,rr

      save
! ----------------------------------------------------------------------
      ! TC-10: t1900 read via state%timecontrol tc_t1900 alias.
      ! [SS-BMI2 Task 4] tstart added to associate
      ! [SS-GR-ATM B.5] at_tav alias for tav read migration
      associate( tc_t1900 => state%timecontrol%t1900, &  ! TC-10
                 tstart   => state%timecontrol%tstart, &  ! [SS-BMI2 Task 4]
                 at_tav   => state%atmosphere%Tav      )  ! [SS-GR-ATM B.5]

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
               if (state%crop%common%icrop >= 1 .and. state%crop%common%icrop <= size(crop_config_global%rotation_loaded)) then  ! [GR-CROPWS B1]
                  if (crop_config_global%rotation_loaded(state%crop%common%icrop)) use_cache = .true.  ! [GR-CROPWS B1]
               end if
            end if
         end if
         if (use_cache) then
            call cropfixed_init_from_config(crop_config_global%rotation_fixed(state%crop%common%icrop), state%crop%common%icrop, lcc, state)  ! [GR-CROPWS B1]
            ! [SS-GR-ATM A5.1] state passed for dual-write of kdif/kdir/swcf/cofab
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
      if (state%crop%common%swrd.eq.1) then
        state%crop%common%rdm = rdmax
      else
        state%crop%common%rdm = min(rdmax,state%crop%common%rdc)   ! [GR-CROPWS B1] state%crop%common%rdc → state%crop%common%rdc
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (tc_t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.           &
     &  dabs(tc_t1900 - state%crop%common%cropstart) .lt. tiny) then   ! [GR-CROPWS B1] cropstart(icrop) → state%crop%common%cropstart

        state%crop%common%dvs = 0.0d0

! --- actual rooting depth
        if (state%crop%common%swrd.eq.1) then
          state%crop%common%rd = afgen (rdtb,22,state%crop%common%dvs)                    ! dvs is 0.0 here (just assigned), no separate read needed
          state%crop%common%rd = min(state%crop%common%rd,state%crop%common%rdm)                            ! [GR-CROPWS B1] state%crop%common%rdm → state%crop%common%rdm
        else
          state%crop%common%rd = min(state%crop%common%rdi,state%crop%common%rdm)         ! [GR-CROPWS B1] state%crop%common%rdi, state%crop%common%rdm → state%crop%common%X
        endif
        state%crop%common%rdpot = state%crop%common%rd

      endif

! --- initial lai or sc
      state%crop%lai = afgen (gctb,(2*magrs),state%crop%common%dvs)               ! [GR-CROPWS B1] dvs → state%crop%common%dvs
      if (swgc.eq.2) then
        gc  = state%crop%lai
        state%crop%lai = state%crop%lai*3.0d0
      endif

! --- initial crop factor or crop height
      state%crop%common%cf = afgen (state%crop%fixed%cftb,(2*magrs),state%crop%common%dvs)    ! [GR-CROPWS B1]
      state%crop%common%ch = afgen (state%crop%fixed%chtb,(2*magrs),state%crop%common%dvs)    ! [GR-CROPWS B1]
      if (state%crop%swcf.eq.3) then                                        ! [GR-CROPWS B1] swcf → state%crop%swcf
        cfeic = afgen (state%crop%fixed%cfeictb,(2*magrs),state%crop%common%dvs)  ! [GR-CROPWS B1]
      endif
      if (state%crop%swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1] [GR-CROPWS B1]

! --- initial storage on canopy
      if (swinter.eq.3) then
        state%atmosphere%siccapact = siccaplai*state%crop%lai                                      ! state%crop%lai local (just computed above)
      endif

! --- initial dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),state%crop%common%dvs)        ! [GR-CROPWS B1]

! --- initial ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),state%crop%common%dvs) ! [GR-CROPWS B1]

! --- initialize matric flux potential (SS-CRP C-2.5: hroot/hleaf/mfluxtable
!     init moved to CropGrowth dispatcher which has access to state).
      if (swdrought .eq. 2) then
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,state%mesh%numnod  ! [GR-BH C7]
         twilt(i) = watcon(wiltpoint, &
                            state%soilwater%vg_params(i), &
                            state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                            i, state%soilwater)                    ! [SS-GR-UTILS Task 5]
        enddo
      endif

      return

      case (2)
      continue

! === calculate potential rate and state variables ======================
      case (3)

! === calculate actual rate and state variables ======================

! --- increase in temperature sum
      dtsum = max (0.0d0,at_tav-tbase)  ! [SS-GR-ATM B.5] tav→state%atmosphere%Tav

! --- development rate
      if (idev.eq.1) then
        dvr = 2.0/lcc
      elseif (idev.eq.2) then
        if (state%crop%common%dvs.lt.1.0d0) then                       ! [GR-CROPWS B1] dvs → state%crop%common%dvs
          dvr = dtsum/tsumea
        else
          dvr = dtsum/tsumam
        endif
      endif

! --- water stress
      ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
      if(dabs(state%atmosphere%ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(min(state%soilwater%tra/state%atmosphere%ptra,1.0d0),0.0d0)  ! [SS-SWC S-2.7]
      endif

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      state%crop%common%dvs = min(state%crop%common%dvs+dvr,2.d0)      ! [GR-CROPWS B1] RHS dvs → state%crop%common%dvs
      state%crop%common%tsum = state%crop%common%tsum + dtsum          ! [GR-CROPWS B1]

! --- leaf area index or soil cover fraction
      state%crop%lai = afgen (gctb,(2*magrs),state%crop%common%dvs)               ! [GR-CROPWS B1]
      if (swgc.eq.2) then
        gc  = state%crop%lai
        state%crop%lai = state%crop%lai*3.0d0
      endif

! --- crop factor or crop height
      state%crop%common%cf        = afgen (state%crop%fixed%cftb,(2*magrs),state%crop%common%dvs)    ! [GR-CROPWS B1]
      state%crop%common%ch        = afgen (state%crop%fixed%chtb,(2*magrs),state%crop%common%dvs)    ! [GR-CROPWS B1]
      if (state%crop%swcf.eq.3) then                                              ! [GR-CROPWS B1]
        cfeic     = afgen (state%crop%fixed%cfeictb,(2*magrs),state%crop%common%dvs)  ! [GR-CROPWS B1]
      endif
      if (state%crop%swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1] [GR-CROPWS B1]

! --- update canopy storage capacity
      if (swinter.eq.3) then
        state%atmosphere%siccapact = siccaplai*state%crop%lai                                      ! state%crop%lai local (just computed above)
      endif

! --- dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),state%crop%common%dvs)        ! [GR-CROPWS B1]

! --- ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),state%crop%common%dvs) ! [GR-CROPWS B1]

      case (4)
          
! --- root extension
      if (state%crop%common%swrd.eq.1) then
        state%crop%common%rdpot = afgen (rdtb,22,state%crop%common%dvs)                  ! [GR-CROPWS B1]
        state%crop%common%rdpot = min(state%crop%common%rdpot,state%crop%common%rdm)                       ! [GR-CROPWS B1]
        state%crop%common%rd    = state%crop%common%rdpot
      else
        rrpot = min (state%crop%common%rdm-state%crop%common%rdpot,state%crop%common%rri)  ! [GR-CROPWS B1]
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if (state%atmosphere%ptra.lt.nihil) rrpot = 0.0d0
        state%crop%common%rdpot = state%crop%common%rdpot + rrpot                        ! [GR-CROPWS B1] RHS state%crop%common%rdpot → state%crop%common%rdpot

        rr = min (state%crop%common%rdm-state%crop%common%rd,state%crop%common%rri)  ! [GR-CROPWS B1]
        if (state%atmosphere%ptra.lt.nihil .or.             &
     &      state%soilwater%flWrtNonox) rr = 0.0d0   ! [SS-GR-CROPWS A2] present(state) guard removed
        if (state%crop%common%swdmi2rd.eq.1 .and. state%atmosphere%ptra.ge.nihil) rr = rr * state%soilwater%tra/state%atmosphere%ptra  ! [SS-SWC S-2.7]
        state%crop%common%rd = state%crop%common%rd + rr                                 ! [GR-CROPWS B1]
      endif

      return

      case default
         call fatalerr_collected ('CropFixed', 'Illegal value for TASK')
      end select

      end associate  ! tc_t1900 => state%timecontrol [TC-10]
      return
      end subroutine cropfixed

      end module cropfixed_runtime_mod
