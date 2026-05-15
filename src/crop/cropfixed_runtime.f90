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
! [SS-GR-CROPRT B2] DEFERRED — cropfixed: all remaining variables globals:
!   magrs: array dim (could → swap_array_dimensions, deferred with rest)
!   icrop, dvs, idev, lai, tsum, cf, ch, rd, rdpot: runtime state (dual-write to state%crop%;
!     global still needed pending Phase C global retirement)
!   max_resp_factor: config field (config%crop%fixed%), needs config threading
!   swrd, rdi, rri, rdc, swgc, swcf, swinter, swdrought, swdmi2rd: switches, no state home
!   cropstart, tbase, tsumea, tsumam, rdmax, rdm: config params, no state home
!   siccapact, siccaplai: state%atmosphere%siccapact migrated; siccaplai no home
!   w_root_ss, wiltpoint, twilt, flhydrlift: JvL params, no state home
!   cftb/chtb/cfeictb: state%crop%fixed homes exist (A5.2 dual-write) but migration to
!     state reads deferred — cropfixed takes optional state, substitution needs
!     non-optional refactor or present() guards; deferred to Phase C cleanup
!   gc, cfeic, gctb, rdtb, mrftb, wrtb: fixed-crop tables/scalars, no state home
!   swinco, reltr: switches, no state home
! ----------------------------------------------------------------------
      use variables, only: magrs, icrop, dvs, idev, lai, tsum, cf, ch, &  ! [SS-GR-CROPRT B2] DEFERRED
                           rd, rdpot, max_resp_factor, swrd, rdi, rri,  &
                           rdc, swgc, swcf, swinter, swdrought, swdmi2rd, &
                           cropstart, tbase, tsumea, tsumam, rdmax, rdm, &
                           siccapact, siccaplai, w_root_ss, wiltpoint,   &
                           twilt, flhydrlift, gc, cfeic,                 &
                           gctb, cftb, chtb, cfeictb, rdtb, mrftb, wrtb, &
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
               if (icrop >= 1 .and. icrop <= size(crop_config_global%rotation_loaded)) then
                  if (crop_config_global%rotation_loaded(icrop)) use_cache = .true.
               end if
            end if
         end if
         if (use_cache) then
            call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop, lcc, state)
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
      if (swrd.eq.1) then
        rdm = rdmax
      else
        rdm = min(rdmax,rdc)
      endif
      state%crop%common%rdm = rdm   ! [SS-GR-CROP A5.1]

! --- skip next initialization if crop parameters are read from *.END file
      if (tc_t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.           &
     &  dabs(tc_t1900 - cropstart(icrop)) .lt. tiny) then

        dvs = 0.0d0

! --- actual rooting depth
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,dvs)
          rd = min(rd,rdm)
        else
          rd = min(rdi,rdm)
        endif
        rdpot = rd
        state%crop%common%dvs   = dvs    ! [SS-GR-CROP A5.1]
        state%crop%common%rd    = rd     ! [SS-GR-CROP A5.1]
        state%crop%common%rdpot = rdpot  ! [SS-GR-CROP A5.1]

      endif

! --- initial lai or sc
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif
      state%crop%lai = lai   ! [SS-GR-ATM A5.2] dual-write

! --- initial crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic = afgen (cfeictb,(2*magrs),dvs)
      endif
      state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
      state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
      if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]

! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
        state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
      endif

! --- initial dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- initial ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

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
        if (dvs.lt.1.0d0) then
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
      dvs = min(dvs+dvr,2.d0)
      tsum = tsum + dtsum
      state%crop%common%dvs  = dvs    ! [SS-GR-CROP A5.1]
      state%crop%common%tsum = tsum   ! [SS-GR-CROP A5.1]

! --- leaf area index or soil cover fraction
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif
      state%crop%lai = lai   ! [SS-GR-ATM A5.2] dual-write

! --- crop factor or crop height
      cf        = afgen (cftb,(2*magrs),dvs)
      ch        = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic     = afgen (cfeictb,(2*magrs),dvs)
      endif
      state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
      state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
      if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]

! --- update canopy storage capacity
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
        state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
      endif

! --- dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

      case (4)
          
! --- root extension
      if (swrd.eq.1) then
        rdpot = afgen (rdtb,22,dvs)
        rdpot = min(rdpot,rdm)
        rd    = rdpot
      else
        rrpot = min (rdm-rdpot,rri)
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if (state%atmosphere%ptra.lt.nihil) rrpot = 0.0d0
        rdpot = rdpot + rrpot

        rr = min (rdm-rd,rri)
        if (state%atmosphere%ptra.lt.nihil .or.             &
     &      state%soilwater%flWrtNonox) rr = 0.0d0   ! [SS-GR-CROPWS A2] present(state) guard removed
        if (swdmi2rd.eq.1 .and. state%atmosphere%ptra.ge.nihil) rr = rr * state%soilwater%tra/state%atmosphere%ptra  ! [SS-SWC S-2.7]
        rd = rd + rr
      endif
      state%crop%common%rdpot = rdpot   ! [SS-GR-CROP A5.1]
      state%crop%common%rd    = rd      ! [SS-GR-CROP A5.1]

      return

      case default
         call fatalerr_collected ('CropFixed', 'Illegal value for TASK')
      end select

      end associate  ! tc_t1900 => state%timecontrol [TC-10]
      return
      end subroutine cropfixed

      end module cropfixed_runtime_mod
