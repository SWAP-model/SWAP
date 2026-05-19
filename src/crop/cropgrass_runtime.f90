! cropgrass_runtime.f90
! GR-CROPWS Phase 0 Commit 0.3: grass extracted from cropgrowth.f90.
! Pure relocation — no behavior change.
! [GR-CROPWS B4]: icrop reads → state%crop%common%icrop (task-1 block, 5 sites);
!   cropstart(icrop) at swinco=3 skip check → state%crop%common%cropstart.
!   icrop and cropstart removed from use variables.
!   daycrop NOT migrated (InitializeCrop zeroes global, state not mirrored there).
! ----------------------------------------------------------------------
      module cropgrass_runtime_mod
      implicit none
      private

      public :: grass

      contains

! ----------------------------------------------------------------------
      subroutine grass(task, tsoil, state)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : detailed grass growth routine
! SS-HEAT pre-Task-8: tsoil(:) non-optional dummy arg; callers pass
!   state%heat%tsoil. Threaded through to sumttd calls.
! SS-CRP C-2.5: state added (optional, intent in) to read flWrtNonox.
! SS-TC TC-10: t1900,daynr read via state%timecontrol tc_* aliases;
!   state threaded to sumttd for its own TC reads.
! SS-GR-ATM A5.1: intent changed inout to allow dual-write in cropgrass_init_from_config.
! [SS-GR-CROPWS A4]: state optional removed — all callers pass state; all if(present(state)) guards dropped.
! [GR-CROP Phase B/7] narrow use variables
! ----------------------------------------------------------------------
      ! [SS-GR-CROPRT B8] DEFERRED — grass: all remaining variables globals:
      !   magrs, macp: array dims (could → swap_array_dimensions, deferred with rest)
      !   icrop, dvs, tsum, daycrop: computed in grass loop; dual-write to state%crop%common%
      !     but global canonical pending Phase C
      !   rid, tbase, tdwi, swinco: config params, no state home
      !   wlv/wst/wrt/wso pools + dwlv/dwst/dwrt: computed in grass loop (WOFOST biomass)
      !   wrtmax, wrtmin: root biomass bounds, no state home
      !   cf, ch, cfeic, lai, laipot, laiem, laiexp/pot, laimax: computed in grass loop
      !   cftb/chtb/cfeictb: state%crop%fixed homes (A5.2 dual-write) but grass has optional
      !     state — same constraint as wofost B7 / cropfixed B2; deferred to Phase C
      !   rdtb, slatb, rgrlai etc.: no state home; rd/rdpot etc.: computed in loop
      !   config switches (swrd etc.), physiology params (reltr, cvl etc.): no state home
      !   leaf arrays (lv/lvpot etc.), JvL params (twilt etc.): no state home
      !   cropstartact/endact/pot: state%crop%grass homes (A5) but written here — Phase C
      !   cuptgraz/pot, tagp/pot, tagpt/pot, seqgrazmow/pot, mowrest, dateharvest:
      !     state%crop%grass/common homes (A5) but written in grass loop — Phase C
      !   pgass/pgasspot: state%crop%wofost homes (A4 dual-write) but written here — Phase C
      !   perdl, dateharvest, lsda: output + harvest tracking, no state home
      !   tsoil: config-staging buffer, renamed to avoid clash with dummy arg
      use variables, only: &                                            ! [SS-GR-CROPRT B8] [GR-CROPWS B4]
        magrs, macp, rid, tbase, daycrop, tdwi, swinco, &  ! [GR-CROPWS B4] icrop/dvs/state%crop%common%tsum retired
        wrtmax, wrtmin,           &  ! wst/wlv/wrt retired
        dwlv, dwlvpot, dwrt, dwrtpot, dwst, dwstpot,                     &
        cf, ch, cfeic, laiem, laiexp, laiexppot, laimax,    &  ! lai/laipot retired
        cftb, chtb, cfeictb, rdtb, slatb, rgrlai, rlwtb, rfsetb,        &
        frtb, fltb, fstb, rdrrtb, rdrstb, kdif,                         &
        rdm, rdmax, rdi, rri, rdc, swrd, swrdc, swdmi2rd,    &  ! rd/rdpot retired
        swdrought, swcf, swgc, swinter, reltr,                           &
        cvl, cvr, cvs, q10, rmr, rml, rms, span, ssa, glaiex, glaiexpot, &
        lv, lvpot, lvage, lvagepot, sla, slapot, ilvold, ilvoldpot,     &
        twilt, wiltpoint, gwrt, siccapact, siccaplai,                   &
        cropstartact, cropendact, cropstartpot, cropendpot,             &
        idaysgraz, idaysgrazpot, idregr, idregrpot,                      &  ! [GR-CROPWS B4] cropstart removed (→state%crop%common%cropstart)
        flgrazing, flgrazingpot, flharvest, flharvestpot,               &
        flhrvendact, flhrvendpot, flhydrlift,                           &
        daygrowth, daygrowthpot, grzdm, dewrest,                        &
        cuptgraz, cuptgrazpot, tagpt, tagptpot,          &  ! tagp retired (→state%crop%wofost%tagp)
        seqgrazmow, seqgrazmowpot, swtsum, iseqgm, iseqgmpot,           &
        iharvest, dmgrztb, dmmowtb, daysgrazingtab, uptgrazingtab,      &
        lossgrazingtab, lossgrztab, lossmowtab,                         &
        delayregrowthtab, zgrz, zmow,                                   &
        mowdm, mowrest, lossdm, plossdm, pmowdm, pgrzdm, pgass, &
        perdl, dateharvest, lsda,                                        &
        dummy_tsoil_gr_ => tsoil
      !! Rename config-staging tsoil to avoid clash with dummy arg tsoil.
      !! [SS-HEAT] Task 9: tsoil retained as config-staging buffer; global is not compute state.
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      ! GR-CROPWS Phase 0: sumttd, update_rootdistribution extracted to cropgrowth_helpers_mod
      use cropgrowth_helpers_mod, only: sumttd, update_rootdistribution

      implicit none

      type(swap_state_t), intent(inout) :: state   ! [SS-GR-CROPWS A4] removed optional — all callers pass state

      integer   i1,task
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil, passed from CropGrowth.
      integer   idelaypot,idelay,i,swhydrlift

      real(8)   laicr,lasum,mres,grazlivinglv,grazlivinglvpot
      real(8)   admi,asrc,ccheck,cvf
      real(8)   dalv,delt,dmi
      real(8)   drrt,drst,dslv,dslv1,dslv2,dslvt,dteff
      real(8)   fcheck,fl,fr,fs,drlv
      real(8)   fysdel,gass,gla,glasol,grlv,grrt,grst
      real(8)   gwst,rest,rmres
      real(8)   slat,teff,twlv,twst
      real(8)   lasumpot,drst1,drst2
      real(8)   drst1pot,drst2pot
      real(8)   gasspot,rmrespot,mrespot,asrcpot,dmipot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot,dslv1pot
      real(8)   dslv2pot,dslvpot,restpot,dalvpot,drlvpot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot,tagpspot,tagps
      real(8)   dummy
      real(8)   dmharvest,dmlastharvest,dmgrazing
      real(8)   lsdb(100)
      real(8)   uptgraz,tagprest,lossgraz
      real(8)   uptgrazpot,lossgrazpot
      integer   daylastharvest,swharvest,daysgrazpot,daysgraz
      integer   swdmmow,swdmgrz
      character(len=11) tmp
      character(len=200) messag
      
! --- rooting
      real(8)   rrpot,rr
      
      integer   maxdaymow,maxdaygrz

      logical   flGrassGrowth
      logical   flDewoolingpot,flDewooling
      character(len=11) ::  dateGrassGrowth
      
!     In case of swlossmow = 1 or swlossgrz = 1: check work-ablility
      integer   nodmow,nodgrz
      integer   swlossmow,swlossgrz
      real(8)   fralossmow
      real(8)   fralossgrz
      real(8)   drz1
      logical   flearlyhrvendpot,flearlyhrvendact
      
      parameter (delt=1.0d0)

      save
! ----------------------------------------------------------------------
      ! SS-TC TC-10: t1900,daynr read via state%timecontrol tc_* aliases.
      ! [SS-BMI2 Task 4] tstart added to associate
      ! [SS-GR-ATM B.5] at_tav alias for tav read migration
      associate( &
        tc_t1900 => state%timecontrol%t1900,  &  ! TC-10
        tc_daynr => state%timecontrol%daynr,  &  ! TC-10
        tstart   => state%timecontrol%tstart, &  ! [SS-BMI2 Task 4]
        at_tav   => state%atmosphere%Tav      &  ! [SS-GR-ATM B.5]
      )

      select case (task)
      case (1)

! === initialization at start of crop =========================================

! --- read grass input data: dispatch on per-rotation typed-config cache
!     (ADR 0016). Falls back to legacy reader for rotations whose
!     .crp.toml is not yet authored. Teardown: end of Phase 4 removes
!     the else-branch.
      block
         use crop_config_global_mod, only: crop_config_global
         use cropgrass_init_mod,     only: cropgrass_init_from_config
         logical :: use_cache
         use_cache = .false.
         if (associated(crop_config_global)) then
            if (allocated(crop_config_global%rotation_loaded)) then
               if (state%crop%common%icrop >= 1 .and. state%crop%common%icrop <= size(crop_config_global%rotation_loaded)) then  ! [GR-CROPWS B4] icrop → state%crop%common%icrop
                  if (crop_config_global%rotation_loaded(state%crop%common%icrop)) then  ! [GR-CROPWS B4]
                     ! Defense-in-depth: only dispatch to cache when the schema
                     ! is fully authored (case 4 + case 2 have amaxtb; the
                     ! hupselbrook skeleton does not — Phase 4 will fill it).
                     if (allocated(crop_config_global%rotation_grass(state%crop%common%icrop)%amaxtb)) &  ! [GR-CROPWS B4]
                        use_cache = .true.
                  end if
               end if
            end if
         end if
         if (use_cache) then
            associate(cfg => crop_config_global%rotation_grass(state%crop%common%icrop))  ! [GR-CROPWS B4] icrop → state%crop%common%icrop
               swharvest      = cfg%swharv
               dmharvest      = cfg%dmharvest
               daylastharvest = int(cfg%daylastharvest)
               dmlastharvest  = cfg%dmlastharvest
               swdmmow        = cfg%swdmmow
               maxdaymow      = cfg%maxdaymow
               swlossmow      = cfg%swlossmow
               swlossgrz      = cfg%swlossgrz
               swdmgrz        = cfg%swdmgrz
               maxdaygrz      = cfg%maxdaygrz
               dmgrazing      = cfg%dmgrazing
               LSDb           = 0.0d0   ! grazing stub-guarded; populated via daysgrazingtab/uptgrazingtab/lossgrazingtab by init
               tagprest       = cfg%tagprest
               swhydrlift     = 0       ! swdrought=2 stub-errored; mirror cropfixed/cropwofost default
               call cropgrass_init_from_config(cfg, state%crop%common%icrop, &  ! [GR-CROPWS B4] icrop → state%crop%common%icrop
                  state%timecontrol%tend, state%timecontrol%tstart, state)  ! [SS-BMI2 Task 4] [SS-GR-ATM A5.1]
            end associate
         else
            ! ADR 0016 cache-miss: typed config required for type=3 rotations.
            ! No silent legacy fallback — the user must author cropgrass.crp.toml.
            call fatalerr_collected('cropgrowth/Grass', &
               'cropgrass rotation has no loaded .crp.toml — author the file or use the legacy executable.')
         end if
      end block

! --- sequence of harvest by mowing, dewooling and grazing
      seqgrazmowpot = seqgrazmow
      state%crop%grass%seqgrazmowpot = seqgrazmowpot   ! [SS-GR-CROP A5.1]

! --- development stage (not used by Grassland, instead Daynrs are used)
      state%crop%common%dvs = -99.99d0

! --- maximum rooting depth
      if (swrd.eq.1) then
        rdm = rdmax
      elseif (swrd.eq.2) then
        rdm = min(rdmax,rdc)
      elseif (swrd.eq.3) then
        rdc = afgen (rlwtb,22,wrtmax)
        rdm = min(rdmax,rdc)
      endif
      state%crop%common%rdm = rdm   ! [SS-GR-CROP A5.1]

! --- skip next initialization if crop parameters are read from *.END file
      if (tc_t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.          &
     &   dabs(tc_t1900 - state%crop%common%cropstart) .lt. tiny) then   ! [GR-CROPWS B4] cropstart(icrop) → state%crop%common%cropstart

        iseqgm = 1
        iseqgmpot = iseqgm

! ---   initial values of crop parameters
        rid = dble(daycrop)
        fr = afgen (frtb,30,rid)
        fl = afgen (fltb,30,rid)
        fs = afgen (fstb,30,rid)
        sla(1) = afgen (slatb,30,rid)
        lvage(1) = 0.d0
        ilvold = 1
        idregr = 0
        slapot(1) = afgen (slatb,30,rid)
        lvagepot(1) = 0.d0
        ilvoldpot = 1
        idregrpot = 0

! ---   initial state variables of the crop
        state%crop%wofost%wrt = fr*tdwi
        wrtmin = state%crop%wofost%wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        state%crop%wofost%wrtpot = state%crop%wofost%wrt
        state%crop%wofost%wst = fs*(1.0d0-fr)*tdwi
        state%crop%wofost%wstpot = state%crop%wofost%wst
        state%crop%wofost%wlv = laiem/sla(1)
        state%crop%wofost%wlvpot = state%crop%wofost%wlv
        
!     KRO-BOO-20160403: intro because comparison with Wofost
        laiem = state%crop%wofost%wlv*sla(1)  ! is not input !
        lv(1) = state%crop%wofost%wlv
        lvpot(1) = lv(1)
        lasum = laiem
        lasumpot = lasum     
        glaiex = 0.0d0
        glaiexpot = 0.0d0
        laiexp = laiem
        laiexppot = laiem
        laimax = laiem
        state%crop%lai = lasum+ssa*state%crop%wofost%wst
        state%crop%common%laipot = state%crop%lai
        dwrt = 0.d0
        dwrtpot = dwrt
        dwlv = 0.d0
        dwlvpot = dwlv
        dwst = 0.d0
        dwstpot = dwst

        daygrowth    = 0
        daygrowthpot = 0

! ---   actual rooting depth
        if (swrd.eq.1) then
          state%crop%common%rd = afgen (rdtb,22,rid)
          state%crop%common%rd = min(state%crop%common%rd,rdm)
        elseif (swrd.eq.2) then
          state%crop%common%rd = min(rdi,rdm)
        elseif (swrd.eq.3) then
          rdi = afgen (rlwtb,22,state%crop%wofost%wrt)
          state%crop%common%rd = min(rdi,rdm)
        endif
        state%crop%common%rdpot = state%crop%common%rd
        
! ---   initial summation variables of the crop
        state%crop%wofost%tagp = state%crop%wofost%wlv+state%crop%wofost%wst
        state%crop%wofost%tagppot = state%crop%wofost%tagp
        tagpt = 0.0d0
        tagptpot = 0.0d0
        cuptgraz = 0.0d0
        cuptgrazpot = 0.0d0
        state%crop%common%tsum = 0.0d0
        
        cropstartpot     = rid
        cropstartact     = rid
        flhrvendpot      = .false.
        flearlyhrvendpot = .false.
        ! [SS-GR-CROP A5.1] mirror grass init-time state
        state%crop%wofost%dwrt        = dwrt
        state%crop%wofost%dwrtpot     = dwrtpot
        state%crop%wofost%dwlv        = dwlv
        state%crop%wofost%dwlvpot     = dwlvpot
        state%crop%wofost%dwst        = dwst
        state%crop%wofost%dwstpot     = dwstpot
        state%crop%wofost%tagpt       = tagpt
        state%crop%wofost%tagptpot    = tagptpot
        state%crop%common%cuptgraz    = cuptgraz
        state%crop%common%cuptgrazpot = cuptgrazpot
        state%crop%grass%cropstartpot = cropstartpot
        state%crop%grass%cropstartact = cropstartact
        
        if (swtsum.eq.0) then
          flGrassGrowth = .true.
        else
          flGrassGrowth = .false.  
        endif
        if (swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          call sumttd('initial',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif

! --- end skip above initialization if crop parameters are read from *.END file
      endif

      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),rid)
        ch = afgen (chtb,(2*magrs),rid)
      else
        cf        = afgen (cftb,(2*magrs),state%crop%lai)
        cfeic     = afgen (cfeictb,(2*magrs),state%crop%lai)
        ch        = afgen(chtb,(2*magrs),state%crop%lai)
      endif
      state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
      state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
      if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]

! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*state%crop%lai
        state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
      endif

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

! --- harvest
!     initialise 
      if (swharvest.eq.2) then
        iharvest = 1
        do while (tc_t1900 .gt. dateharvest(iharvest))
          iharvest = iharvest + 1
        enddo
      endif      
      
! --- Find node for monitoring work-ability
      if (swlossmow .eq. 1) then
         
        ! Find node for monitoring work-ability during mowing
        nodmow = 1
        drz1       = -1.d0 * zmow - state%mesh%dz(nodmow)  ! [GR-BH C7]
        do while (drz1 .gt. 0.d0)
          nodmow   = nodmow + 1
          drz1 = drz1 - state%mesh%dz(nodmow)  ! [GR-BH C7]
        enddo
      
      endif   
      
      if (swlossgrz .eq. 1) then
        
        ! Find node and layer for monitoring work-ability at start of grazing
        nodgrz = 1
        drz1       = -1.d0 * zgrz - state%mesh%dz(nodgrz)  ! [GR-BH C7]
        do while (drz1 .gt. 0.d0)
          nodgrz   = nodgrz + 1
          drz1 = drz1 - state%mesh%dz(nodgrz)  ! [GR-BH C7]
        enddo
         
      endif
      
      return

      case (2)

! === calculate potential rate and state variables ======================================

! --- rates of change of the grass variables ---------------------------------------------

      rid = dble(daycrop)
      
! --- check end of harvest
      if (flhrvendpot) then
        if (flearlyhrvendpot) then
          cropstartpot  = rid - 1.d0
        else
          cropstartpot  = rid
        endif
        pmowdm        = 0.d0
        pgrzdm        = 0.d0
        plossdm       = 0.d0
        state%crop%grass%cropstartpot = cropstartpot   ! [SS-GR-CROP A5.1]
        state%crop%wofost%plossdm     = plossdm        ! [SS-GR-CROP A5.1]
      endif
      flhrvendpot      = .false.
      flearlyhrvendpot = .false.

! --- grass growth initiated by tsum from 1st day of calendar year
      state%crop%common%tsum = state%crop%common%tsum + max(0.0d0,at_tav)  ! [SS-GR-ATM B.5]
      if (.not. flGrassGrowth) then
        
        ! grass growth initiated by tsum
        if (swtsum.eq.1) then
          if (state%crop%common%tsum.ge.200.d0) then
            flGrassGrowth = .true.
          endif
        endif
        
        ! grass growth initiated by temperature, time and depth
        if (swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          if (dateGrassGrowth.eq.'undefined') call sumttd('dynamic',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif
      
        ! check if grass growth has started
        if (flGrassGrowth) then
          cropstartpot = rid
          cropstartact = rid
        endif

      endif
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.idregrpot) then

! ===   daily dry matter production ===

        gasspot = state%crop%wofost%pgasspot

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(rmr*state%crop%wofost%wrtpot+rml*state%crop%wofost%wlvpot+rms*state%crop%wofost%wstpot)*afgen(rfsetb,30,rid)
        teff = q10**((at_tav-25.0d0)/10.0d0)  ! [SS-GR-ATM B.5]
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)
! ---   check on partitioning
        fcheck = fr+(fl+fs)*(1.0d0-fr) - 1.0d0
        if (dabs(fcheck).gt.0.0001d0) then
          write(tmp,'(f6.3)') rid
          tmp = adjustl (tmp)
          Messag ='The sum of partitioning factors for leaves, stems'// &
     &    ' and storage organs is not equal to one at time '            &
     &    //trim(tmp)//'.'
          call fatalerr_collected ('grass_pot',messag)
        endif

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmipot = cvf*asrcpot

! ---   check on carbon balance
        ccheck = (gasspot-mrespot-(fr+(fl+fs)*(1.0d0-fr))*dmipot/cvf)   &
     &         /max(0.0001d0,gasspot)
        if (dabs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr_collected ('grass_pot',messag)
        endif


! ===   growth rate by plant organ ===

! ---   growth rate roots and aerial parts

        grrtpot = fr*dmipot
        ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
        ! growth of the roots is balanced by the death of root tissue
        if (swrd.eq.3 .and. state%crop%wofost%wrtpot.gt.wrtmax) then
          drrtpot = grrtpot
          drrtpot = max(drrtpot,state%crop%wofost%wrtpot*afgen (rdrrtb,30,rid))
        else  
          drrtpot = state%crop%wofost%wrtpot*afgen (rdrrtb,30,rid)
        endif  
        gwrtpot = grrtpot - drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/kdif
        dslv2pot=state%crop%wofost%wlvpot*max(0.0d0,                                      &
     &                  min(0.03d0,0.03d0*(state%crop%common%laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = ilvoldpot

        do while (restpot.gt.lvpot(max(i1,1)).and.i1.ge.1)
          restpot = restpot-lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (lvagepot(max(i1,1)).gt.span.and.restpot.gt.0.and.           &
     &                          i1.ge.1) then
          dalvpot = lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvagepot(max(i1,1)).gt.span)
          dalvpot = dalvpot+lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot   = dslvpot+dalvpot

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (slatb,30,rid)
        if (laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,at_tav-tbase)  ! [SS-GR-ATM B.5]
          glaiexpot = laiexppot*rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (rdrstb,30,rid)*state%crop%wofost%wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        daygrowthpot = daygrowthpot + 1

!       Check trigger to start mowing event
        if (seqgrazmowpot(iseqgmpot) .eq. 2) then

          flharvestpot = .false.
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (state%crop%wofost%tagppot .gt. dmharvest .or. (tc_daynr .gt. daylastharvest  &
     &          .and. state%crop%wofost%tagppot .gt. dmlastharvest)) then
                flharvestpot = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(dmmowtb,20,rid)
              if (state%crop%wofost%tagppot .gt. dmharvest .or.                           &
     &                 (daygrowthpot .gt. maxdaymow .and. iseqgmpot .gt. 1)) then
                flharvestpot = .true.
              endif
            endif

          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(tc_t1900 .gt. dateharvest(iharvest)) then
              flharvestpot = .true.
            endif
          endif

!         In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
          if (flharvestpot) then
            iseqgmpot = iseqgmpot + 1
            slapot(1) = afgen (slatb,30,rid)
            fl = afgen (fltb,30,rid)
            fs = afgen (fstb,30,rid)
            state%crop%wofost%wlvpot = mowrest / (1.d0 + (fs/fl))
            state%crop%wofost%wstpot = fs/fl*state%crop%wofost%wlvpot
            dwlvpot = 0.0d0
            dwstpot = 0.0d0
            lvagepot(1) = 0.0d0
            ilvoldpot = 1
            lasumpot = state%crop%wofost%wlvpot * slapot(1)
            laiexppot = lasumpot
            lvpot(1) = state%crop%wofost%wlvpot
            
            gwstpot = 0.0d0
            gwrtpot = 0.0d0
            drlvpot = 0.0d0
            drstpot = 0.0d0
            drrtpot = 0.0d0
            
            daygrowthpot = 0
            
!           losses due to treading
            fralossmow = 0.d0
            if (swlossmow.eq.1) then
              fralossmow = afgen(lossmowtab,200,state%soilwater%h(nodmow))  ! [SS-SWC S-2.7]
            end if

!           harvest
            tagpspot = max(0.0d0,(state%crop%wofost%tagppot-(state%crop%wofost%wlvpot+dwlvpot+state%crop%wofost%wstpot+dwstpot)))
            tagptpot = tagptpot + tagpspot * (1.d0 - FraLossMow)

            cropendpot  = rid
            flhrvendpot = .true.
            pmowdm      = tagpspot * (1.d0 - FraLossMow)
            plossdm     = tagpspot * FraLossMow
            
!           set regrowth delay
            idelaypot = int(afgen(DelayRegrowthTab,200,tagpspot))
            idregrpot = daycrop + idelaypot

          endif          
          
!       Check trigger to start grazing event          
        else if (seqgrazmowpot(iseqgmpot) .eq. 1 .or. seqgrazmowpot(iseqgmpot) .eq. 3) then

          flharvestpot = .false.
          
          if (.not. flgrazingpot) then
              
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (state%crop%wofost%tagppot .gt. dmgrazing) then
                  flharvestpot = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(dmgrztb,20,rid)
                if (state%crop%wofost%tagppot .gt. dmgrazing .or.                           &
     &                 (daygrowthpot .gt. maxdaygrz .and. iseqgmpot .gt. 1)) then
                  flharvestpot = .true.
                endif
              endif
              
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(tc_t1900 .gt. dateharvest(iharvest)) then
                flharvestpot = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (flharvestpot .or. flgrazingpot) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgrazpot = lsda(iseqgmpot) *                                 &
     &                         afgen(uptgrazingtab,200,lsda(iseqgmpot))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgrazpot = lsda(iseqgmpot) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgmpot))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(lossgrztab,200,state%soilwater%h(nodgrz))  ! [SS-SWC S-2.7]
            end if
            lossgrazpot = lossgrazpot + state%crop%wofost%tagppot * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. flgrazingpot) then
              daygrowthpot   = 0
              idaysgrazpot   = 0
              flDewoolingpot = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((state%crop%wofost%tagppot - uptgrazpot - lossgrazpot) .gt. tagprest) then
              
              flgrazingpot = .true.
              cuptgrazpot  = cuptgrazpot + uptgrazpot
          
!             distribute grazing over stems and leaves (living and dead parts)
              state%crop%wofost%wstpot  = state%crop%wofost%wstpot  - (uptgrazpot+lossgrazpot) * state%crop%wofost%wstpot  / state%crop%wofost%tagppot
              dwstpot = dwstpot - (uptgrazpot+lossgrazpot) * dwstpot / state%crop%wofost%tagppot
              dwlvpot = dwlvpot - (uptgrazpot+lossgrazpot) * dwlvpot / state%crop%wofost%tagppot
              grazlivinglvpot =   (uptgrazpot+lossgrazpot) * state%crop%wofost%wlvpot  / state%crop%wofost%tagppot
          
!             reduce leave weights
              i1 = ilvoldpot
              do while (grazlivinglvpot .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglvpot .ge. lvpot(i1)) then
                  grazlivinglvpot = grazlivinglvpot - lvpot(i1)
                  lvpot(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  lvpot(i1) = lvpot(i1) - grazlivinglvpot
                  grazlivinglvpot = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              cropendpot = rid
              pgrzdm     = pgrzdm + uptgrazpot
              plossdm    = state%crop%wofost%tagppot * fralossgrz
              
!             Check number of days with grazing
              daysgrazpot  = int(afgen(daysgrazingtab,200,lsda(iseqgmpot)))
              idaysgrazpot = idaysgrazpot + 1
              if(idaysgrazpot .eq. daysgrazpot) then
                flgrazingpot = .false.
                flhrvendpot  = .true.
                if (seqgrazmowpot(iseqgmpot) .eq. 3) then
                  flDewoolingpot  = .true.
                endif
                daygrowthpot = 0
                iseqgmpot = iseqgmpot + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazingpot .or. swharvest .eq. 2) then
              flgrazingpot     = .false.
              flhrvendpot      = .true.
              flearlyhrvendpot = .true.
              if (seqgrazmowpot(iseqgmpot) .eq. 3 .and. state%crop%wofost%tagppot .gt. dewrest) then
                flDewoolingpot   = .true.
                flearlyhrvendpot = .false.
              endif
              daygrowthpot = 0
              iseqgmpot = iseqgmpot + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            idregrpot = daycrop

!           Dewooling after grazing event            
            if (flDewoolingpot) then

              slapot(1) = afgen (slatb,30,rid)
              fl = afgen (fltb,30,rid)
              fs = afgen (fstb,30,rid)
              state%crop%wofost%wlvpot = dewrest / (1.d0 + (fs/fl))
              state%crop%wofost%wstpot = fs/fl*state%crop%wofost%wlvpot
              dwlvpot = 0.0d0
              dwstpot = 0.0d0
              lvagepot(1) = 0.0d0
              ilvoldpot = 1
              lasumpot = state%crop%wofost%wlvpot * slapot(1)
              laiexppot = lasumpot
              lvpot(1) = state%crop%wofost%wlvpot
              
              gwstpot = 0.0d0
              gwrtpot = 0.0d0
              drlvpot = 0.0d0
              drstpot = 0.0d0
              drrtpot = 0.0d0
              
!             Assumption: one day delay in regrowth after grazing
              idregrpot = daycrop + 1
              
            endif
            
          endif
          
        endif
        
        if (daycrop .ge. idregrpot) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(at_tav-tbase)/(35.0d0-tbase))  ! [SS-GR-ATM B.5]

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = ilvoldpot
           do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.lvpot(i1)) then
              dslvtpot = dslvtpot-lvpot(i1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              lvpot(i1) = lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvagepot(max(i1,1)) .gt. span .and. i1 .ge. 1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvoldpot,1,-1
            lvpot(i1+1) = lvpot(i1)
            slapot(i1+1) = slapot(i1)
            lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
          enddo
          ilvoldpot = ilvoldpot+1

! ---     new leaves in class 1
          lvpot(1) = grlvpot*delt
          slapot(1) = slatpot
          lvagepot(1) = 0.d0

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          state%crop%wofost%wlvpot = 0.d0
          do i1 = 1,ilvoldpot
            lasumpot = lasumpot+lvpot(i1)*slapot(i1)
            state%crop%wofost%wlvpot = state%crop%wofost%wlvpot+lvpot(i1)
          enddo

          laiexppot = laiexppot+glaiexpot*delt

        endif

! ---   dry weight of living plant organs
        state%crop%wofost%wrtpot = state%crop%wofost%wrtpot+gwrtpot*delt
        state%crop%wofost%wstpot = state%crop%wofost%wstpot+gwstpot*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrtpot = dwrtpot+drrtpot*delt
        dwlvpot = dwlvpot+drlvpot*delt
        dwstpot = dwstpot+drstpot*delt

! ---   dry weight of dead and living plant organs
        twlvpot = state%crop%wofost%wlvpot+dwlvpot
        twstpot = state%crop%wofost%wstpot+dwstpot
        state%crop%wofost%tagppot = twlvpot+twstpot

! ---   leaf area index
        state%crop%common%laipot = lasumpot+ssa*state%crop%wofost%wstpot
!       prevent immediate lai reduction at emergence
!       KRO-BOO-20160403: suppressed because deviates from Wofost
!       laipot = max(laipot, laiem)

        ! root extension
        if (swrd.eq.1) then
          state%crop%common%rdpot = afgen (rdtb,22,rid)
          state%crop%common%rdpot = min(state%crop%common%rdpot,rdm)
        elseif (swrd.eq.2) then
          rrpot = min (rdm-state%crop%common%rdpot,rri)
          if (fr.le.0.0d0 .or. state%crop%wofost%pgasspot.lt.1.0d0) rrpot = 0.0d0
          state%crop%common%rdpot = state%crop%common%rdpot + rrpot
        elseif (swrd.eq.3) then
          state%crop%common%rdpot = afgen (rlwtb,22,state%crop%wofost%wrtpot)
          state%crop%common%rdpot = min(state%crop%common%rdpot,rdm)
        endif

      endif

      ! [SS-GR-CROP A5.1] mirror grass case(2) potential state
      state%crop%wofost%dwrtpot       = dwrtpot
      state%crop%wofost%dwlvpot       = dwlvpot
      state%crop%wofost%dwstpot       = dwstpot
      state%crop%wofost%plossdm       = plossdm
      state%crop%common%cuptgrazpot   = cuptgrazpot
      state%crop%wofost%tagptpot      = tagptpot
      state%crop%grass%cropstartpot   = cropstartpot
      state%crop%grass%cropendpot     = cropendpot

      return

      case (3)

! === calculate actual rate and state variables ======================================

! --- check end of harvest
      if (flhrvendact) then
        if (flearlyhrvendact) then
          cropstartact  = rid - 1.d0
        else
          cropstartact  = rid
        endif
        mowdm        = 0.d0
        grzdm        = 0.d0
        lossdm       = 0.d0
        state%crop%grass%cropstartact = cropstartact   ! [SS-GR-CROP A5.1]
        state%crop%wofost%lossdm      = lossdm         ! [SS-GR-CROP A5.1]
      endif
      flhrvendact      = .false.
      flearlyhrvendact = .false.

! --- rates of change of the crop variables ---------------------------------------------
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop .ge. idregr) then

! ===   daily dry matter production ===

! ---   water stress reduction of pgass to gass
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if(dabs(state%atmosphere%ptra).lt.nihil) then
          reltr = 1.0d0
        else
          reltr = max(0.0d0,min(1.0d0,state%soilwater%tra/state%atmosphere%ptra))  ! [SS-SWC S-2.7]
        endif
        gass = pgass * reltr

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (rmr*state%crop%wofost%wrt+rml*state%crop%wofost%wlv+rms*state%crop%wofost%wst)*afgen(rfsetb,30,rid)
        teff = q10**((at_tav-25.0d0)/10.0d0)  ! [SS-GR-ATM B.5]
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   partitioning factors (relevant for restart)
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmi = cvf*asrc
! ---   check on carbon balance
        ccheck = (gass-mres-(fr+(fl+fs)*(1.0d0-fr))*dmi/cvf)            &
     &         /max(0.0001d0,gass)      
        if (dabs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr_collected ('grass_act',messag)
        endif

! ===   growth rate by plant organ ===

! ---   growth rate roots and aerial parts
        ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
        ! growth of the roots is balanced by the death of root tissue
        grrt = fr*dmi
        if (swrd.eq.3 .and. state%soilwater%flWrtNonox) grrt = 0.d0   ! [SS-GR-CROPWS A4] present(state) guard removed
        if (swrd.eq.3 .and. state%crop%wofost%wrt.gt.wrtmax) then
          drrt = grrt
          drrt = max(drrt,state%crop%wofost%wrt*afgen (rdrrtb,30,rid))
        else  
          drrt = state%crop%wofost%wrt*afgen (rdrrtb,30,rid)
        endif  
        gwrt = grrt-drrt

! ---   growth rate leaves

! ---   weight of new leaves
        admi = (1.0d0-fr)*dmi        
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        dslv1 = state%crop%wofost%wlv*(1.0d0-reltr)*perdl
        laicr = 3.2d0/kdif
        dslv2 = state%crop%wofost%wlv*max(0.0d0,min(0.03d0,0.03d0*(state%crop%lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = ilvold

        do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
          rest = rest-lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (lvage(max(i1,1)).gt.span.and.rest.gt.0.and.i1.ge.1) then
          dalv = lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
          dalv = dalv+lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv   = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (laiexp.lt.6.0d0) then
          dteff = max (0.0d0,at_tav-tbase)  ! [SS-GR-ATM B.5]
          glaiex = laiexp*rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        drst1 = state%crop%wofost%wst*(1.0d0-reltr)*perdl
! ---   death of stems due to ageing
        drst2 = afgen (rdrstb,30,rid)*state%crop%wofost%wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        daygrowth = daygrowth + 1

!       Check trigger to start mowing event
        if (seqgrazmow(iseqgm) .eq. 2) then
          
          flharvest = .false.   
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (state%crop%wofost%tagp .gt. dmharvest .or. (tc_daynr .gt. daylastharvest  &
     &          .and. state%crop%wofost%tagp .gt. dmlastharvest)) then
                flharvest = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(dmmowtb,20,rid)
              if (state%crop%wofost%tagp .gt. dmharvest .or.                           &
     &                 (daygrowth .gt. maxdaymow .and. iseqgm .gt. 1)) then
                flharvest = .true.
              endif
            endif
          
          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(tc_t1900 .gt. dateharvest(iharvest)) then
              iharvest = iharvest + 1
              flharvest = .true.
            endif
          endif
          
!       In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
        if (flharvest) then
          iseqgm = iseqgm + 1
          sla(1) = afgen (slatb,30,rid)
          fl = afgen (fltb,30,rid)
          fs = afgen (fstb,30,rid)
          state%crop%wofost%wlv = mowrest / (1.d0 + (fs/fl))
          state%crop%wofost%wst = fs/fl*state%crop%wofost%wlv
          dwlv = 0.0d0
          dwst = 0.0d0
          lvage(1) = 0.0d0
          ilvold = 1
          lasum = state%crop%wofost%wlv * sla(1)
          laiexp = lasum
          lv(1) = state%crop%wofost%wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          daygrowth = 0

!         losses due to treading
          FraLossMow = 0.d0
          if (swlossmow.eq.1) then
            FraLossMow = afgen(lossmowtab,200,state%soilwater%h(nodmow))  ! [SS-SWC S-2.7]
          end if
          
!         harvest
          tagps = max (0.0d0,(state%crop%wofost%tagp-(state%crop%wofost%wlv+dwlv+state%crop%wofost%wst+dwst)))
          tagpt = tagpt + tagps * (1.d0 - fralossmow)

          cropendact  = rid
          flhrvendact = .true.
          mowdm   = tagps * (1.d0 - FraLossMow)
          lossdm  = tagps * FraLossMow
          
! ---     set regrowth delay
          idelay = int(afgen(DelayRegrowthTab,200,tagps))
          idregr = daycrop + idelay

        endif
          
!       Check trigger to start grazing event          
        else if (seqgrazmow(iseqgm) .eq. 1 .or. seqgrazmow(iseqgm) .eq. 3) then

          flharvest = .false.
            
          if (.not. flgrazing) then
            
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (state%crop%wofost%tagp .gt. dmgrazing) then
                  flharvest = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(dmgrztb,20,rid)
                if (state%crop%wofost%tagp .gt. dmgrazing .or.                           &
     &            (daygrowth .gt. maxdaygrz .and. iseqgm .gt. 1)) then
                  flharvest = .true.
                endif
              endif
            
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(tc_t1900 .gt. dateharvest(iharvest)) then
                iharvest = iharvest + 1
                flharvest = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (flharvest .or. flgrazing) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = lsda(iseqgm)*afgen(uptgrazingtab,200,lsda(iseqgm))            

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = lsda(iseqgm) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgm))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(lossgrztab,200,state%soilwater%h(nodgrz))  ! [SS-SWC S-2.7]
            end if
            lossgraz = lossgraz + state%crop%wofost%tagp * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. flgrazing) then
              daygrowth   = 0
              idaysgraz   = 0
              flDewooling = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((state%crop%wofost%tagp - uptgraz - lossgraz) .gt. tagprest) then
              
              flgrazing = .true.
              cuptgraz  = cuptgraz + uptgraz

!             distribute grazing over stems and leaves (living and dead parts)
              state%crop%wofost%wst  = state%crop%wofost%wst  -  (uptgraz+lossgraz) * state%crop%wofost%wst  / state%crop%wofost%tagp
              dwst = dwst -  (uptgraz+lossgraz) * dwst / state%crop%wofost%tagp
              dwlv = dwlv -  (uptgraz+lossgraz) * dwlv / state%crop%wofost%tagp
              grazlivinglv = (uptgraz+lossgraz) * state%crop%wofost%wlv  / state%crop%wofost%tagp
          
!             reduce leave weights
              i1 = ilvold
              do while (grazlivinglv .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglv .ge. lv(i1)) then
                  grazlivinglv = grazlivinglv - lv(i1)
                  lv(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  lv(i1) = lv(i1) - grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              cropendact = rid
              grzdm      = grzdm + uptgraz
              lossdm     = state%crop%wofost%tagp * fralossgrz
              
!             Check number of days with grazing
              daysgraz  = int(afgen(daysgrazingtab,200,lsda(iseqgm)))
              idaysgraz = idaysgraz + 1
              if(idaysgraz .eq. daysgraz) then
                flgrazing   = .false.
                flhrvendact = .true.
                if (seqgrazmow(iseqgm) .eq. 3) then
                  flDewooling  = .true.
                endif
                daygrowth = 0
                iseqgm = iseqgm + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazing .or. swharvest .eq. 2) then
              flgrazing        = .false.
              flhrvendact      = .true.
              flearlyhrvendact = .true.
              if (seqgrazmow(iseqgm) .eq. 3 .and. state%crop%wofost%tagp .gt. dewrest) then
                flDewooling      = .true.
                flearlyhrvendact = .false.
              endif
              daygrowth = 0
              iseqgm = iseqgm + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            idregr = daycrop

!           Dewooling after grazing event            
            if (flDewooling) then

              sla(1) = afgen (slatb,30,rid)
              fl = afgen (fltb,30,rid)
              fs = afgen (fstb,30,rid)
              state%crop%wofost%wlv = dewrest / (1.d0 + (fs/fl))
              state%crop%wofost%wst = fs/fl*state%crop%wofost%wlv
              dwlv = 0.0d0
              dwst = 0.0d0
              lvage(1) = 0.0d0
              ilvold = 1
              lasum = state%crop%wofost%wlv * sla(1)
              laiexp = lasum
              lv(1) = state%crop%wofost%wlv
    
              gwst = 0.0d0
              gwrt = 0.0d0
              drlv = 0.0d0
              drst = 0.0d0
              drrt = 0.0d0
    
!             Assumption: one day delay in regrowth after grazing
              idregr = daycrop + 1

            endif
            
          endif
          
        endif

        if (daycrop .ge. idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(at_tav-tbase)/(35.0d0-tbase))  ! [SS-GR-ATM B.5]

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.lv(i1)) then
              dslvt = dslvt-lv(i1)
              lv(i1) = 0.0d0
              i1 = i1-1
            else
              lv(i1) = lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvage(max(i1,1)).gt.span.and.i1.ge.1)
              lv(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvold,1,-1
            lv(i1+1) = lv(i1)
            sla(i1+1) = sla(i1)
            lvage(i1+1) = lvage(i1)+fysdel*delt
          enddo
          ilvold = ilvold+1

! ---     new leaves in class 1
          lv(1) = grlv*delt
          sla(1) = slat
          lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          state%crop%wofost%wlv = 0.d0
          do i1 = 1,ilvold
            lasum = lasum+lv(i1)*sla(i1)
            state%crop%wofost%wlv = state%crop%wofost%wlv+lv(i1)
          enddo

          laiexp = laiexp+glaiex*delt

        endif

! ---   dry weight of living plant organs
        state%crop%wofost%wrt = state%crop%wofost%wrt+gwrt*delt
        state%crop%wofost%wst = state%crop%wofost%wst+gwst*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrt = dwrt+drrt*delt
        dwlv = dwlv+drlv*delt
        dwst = dwst+drst*delt

! ---   dry weight of dead and living plant organs
        twlv = state%crop%wofost%wlv+dwlv
        twst = state%crop%wofost%wst+dwst
        state%crop%wofost%tagp = twlv+twst

! ---   leaf area index
        state%crop%lai = lasum+ssa*state%crop%wofost%wst
        laimax = max (state%crop%lai,laimax)

! ---   update normalized cumulative root density based on root extraction or stress (cumdens)
        if (swrdc .eq. 1) call update_rootdistribution(state)
        
        ! root extension
        if (swrd.eq.1) then
          state%crop%common%rd = afgen (rdtb,22,rid)
          state%crop%common%rd = min(state%crop%common%rd,rdm)
        elseif (swrd.eq.2) then
          rr = min (rdm-state%crop%common%rd,rri)
          if (fr.le.0.0d0 .or. pgass.lt.1.0d0 .or.                    &
     &        state%soilwater%flWrtNonox) rr = 0.0d0   ! [SS-GR-CROPWS A4] present(state) guard removed
          if (swdmi2rd.eq.1 .and. pgass.ge.1.0d0)              rr = rr * gass/pgass
          state%crop%common%rd = state%crop%common%rd + rr
        elseif (swrd.eq.3) then
          state%crop%common%rd = afgen (rlwtb,22,state%crop%wofost%wrt)
          state%crop%common%rd = min(state%crop%common%rd,rdm)
        endif

! ---   set crop height and cropfactor
        if (swcf.ne.3) then
          cf = afgen (cftb,(2*magrs),rid)
          ch = afgen (chtb,(2*magrs),rid)
        else
          cf = afgen (cftb,(2*magrs),state%crop%lai)
          cfeic = afgen (cfeictb,(2*magrs),state%crop%lai)
          ch = afgen(chtb,(2*magrs),state%crop%lai)
        endif
        state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
        state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
        if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]

! ---   update canopy storage capacity
        if (swinter.eq.3) then
          siccapact = siccaplai*state%crop%lai
          state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
        endif

      endif

      ! [SS-GR-CROP A5.1] mirror grass case(3) actual state
      state%crop%wofost%dwrt        = dwrt
      state%crop%wofost%dwlv        = dwlv
      state%crop%wofost%dwst        = dwst
      state%crop%wofost%tagp        = state%crop%wofost%tagp
      state%crop%wofost%pgass       = pgass
      state%crop%wofost%lossdm      = lossdm
      state%crop%common%cuptgraz    = cuptgraz
      state%crop%wofost%tagpt       = tagpt
      state%crop%grass%cropstartact = cropstartact
      state%crop%grass%cropendact   = cropendact

      return

      case default
         call fatalerr_collected ('Grass', 'Illegal value for TASK')
      end select

      end associate  ! tc_t1900, tc_daynr => state%timecontrol [TC-10]
      return
      end subroutine grass

      end module cropgrass_runtime_mod
