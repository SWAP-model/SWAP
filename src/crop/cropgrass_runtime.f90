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
      !   leaf arrays (state%crop%common%lv/state%crop%common%lvpot etc.), JvL params (twilt etc.): no state home
      !   cropstartact/endact/pot: state%crop%grass homes (A5) but written here — Phase C
      !   state%crop%common%cuptgraz/pot, tagp/pot, tagpt/pot, seqgrazmow/pot, mowrest, dateharvest:
      !     state%crop%grass/common homes (A5) but written in grass loop — Phase C
      !   pgass/pgasspot: state%crop%wofost homes (A4 dual-write) but written here — Phase C
      !   perdl, dateharvest, state%crop%grass%lsda: output + harvest tracking, no state home
      !   tsoil: config-staging buffer, renamed to avoid clash with dummy arg
      use variables, only: &                                            ! [SS-GR-CROPRT B8] [GR-CROPWS B4]
        magrs, macp, rid, daycrop, swinco,       &  ! tdwi/[GR-CROPWS B4] icrop/dvs/tsum/wst/wlv/wrt/dw*/tbase retired
        wrtmin,                  &  ! wrtmax retired
        ! laiem/laiexp/laiexppot/laimax/lai/laipot/cfeic retired
        ! cftb/chtb/cfeictb/rdtb/rlwtb/slatb/rgrlai/rfsetb retired
        ! frtb/fltb/fstb/rdrrtb/rdrstb retired
        rdmax,                           &  ! rd/rdpot/swrd/swrdc/swdmi2rd retired
        reltr,                                                           &  ! swgc/swdrought/swinter/swcf retired
        ! glaiex/glaiexpot/cvl/cvr/cvs/q10/rmr/rml/rms/span/ssa retired
        ! lv/lvpot/lvage/lvagepot/sla/slapot/ilvold/ilvoldpot retired
        twilt, wiltpoint, gwrt, siccaplai,                   &  ! cropstartact/endact/startpot/endpot retired
        flhydrlift,                                                      &  ! state%crop%grass%idaysgraz*/state%crop%grass%idregr*/state%crop%grass%flGrazing*/state%crop%grass%flHarvest*/flhrvend*/state%crop%grass%daygrowth*/state%crop%grass%grzdm/state%crop%grass%dewrest/state%crop%grass%swtsum/state%crop%grass%iseqgm*/state%crop%grass%iharvest/state%crop%grass%mowdm*/state%crop%grass%pgrzdm/state%crop%grass%pmowdm/state%crop%grass%lsda/state%crop%common%cuptgraz* retired
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
      state%crop%grass%seqgrazmowpot = state%crop%grass%seqgrazmow

! --- development stage (not used by Grassland, instead Daynrs are used)
      state%crop%common%dvs = -99.99d0

! --- maximum rooting depth
      if (state%crop%common%swrd.eq.1) then
        state%crop%common%rdm = rdmax
      elseif (state%crop%common%swrd.eq.2) then
        state%crop%common%rdm = min(rdmax,state%crop%common%rdc)
      elseif (state%crop%common%swrd.eq.3) then
        state%crop%common%rdc = afgen (state%crop%common%rlwtb,22,state%crop%common%wrtmax)
        state%crop%common%rdm = min(rdmax,state%crop%common%rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (tc_t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.          &
     &   dabs(tc_t1900 - state%crop%common%cropstart) .lt. tiny) then   ! [GR-CROPWS B4] cropstart(icrop) → state%crop%common%cropstart

        state%crop%grass%iseqgm = 1
        state%crop%grass%iseqgmpot = state%crop%grass%iseqgm

! ---   initial values of crop parameters
        rid = dble(daycrop)
        fr = afgen (state%crop%common%frtb,30,rid)
        fl = afgen (state%crop%common%fltb,30,rid)
        fs = afgen (state%crop%common%fstb,30,rid)
        state%crop%common%sla(1) = afgen (state%crop%common%slatb,30,rid)
        state%crop%common%lvage(1) = 0.d0
        state%crop%common%ilvold = 1
        state%crop%grass%idregr = 0
        state%crop%common%slapot(1) = afgen (state%crop%common%slatb,30,rid)
        state%crop%common%lvagepot(1) = 0.d0
        state%crop%common%ilvoldpot = 1
        state%crop%grass%idregrpot = 0

! ---   initial state variables of the crop
        state%crop%wofost%wrt = fr*state%crop%common%tdwi
        wrtmin = state%crop%wofost%wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        state%crop%wofost%wrtpot = state%crop%wofost%wrt
        state%crop%wofost%wst = fs*(1.0d0-fr)*state%crop%common%tdwi
        state%crop%wofost%wstpot = state%crop%wofost%wst
        state%crop%wofost%wlv = state%crop%common%laiem/state%crop%common%sla(1)
        state%crop%wofost%wlvpot = state%crop%wofost%wlv
        
!     KRO-BOO-20160403: intro because comparison with Wofost
        state%crop%common%laiem = state%crop%wofost%wlv*state%crop%common%sla(1)  ! is not input !
        state%crop%common%lv(1) = state%crop%wofost%wlv
        state%crop%common%lvpot(1) = state%crop%common%lv(1)
        lasum = state%crop%common%laiem
        lasumpot = lasum     
        state%crop%common%glaiex = 0.0d0
        state%crop%common%glaiexpot = 0.0d0
        state%crop%common%laiexp = state%crop%common%laiem
        state%crop%common%laiexppot = state%crop%common%laiem
        state%crop%common%laimax = state%crop%common%laiem
        state%crop%lai = lasum+state%crop%common%ssa*state%crop%wofost%wst
        state%crop%common%laipot = state%crop%lai
        state%crop%wofost%dwrt = 0.d0
        state%crop%wofost%dwrtpot = state%crop%wofost%dwrt
        state%crop%wofost%dwlv = 0.d0
        state%crop%wofost%dwlvpot = state%crop%wofost%dwlv
        state%crop%wofost%dwst = 0.d0
        state%crop%wofost%dwstpot = state%crop%wofost%dwst

        state%crop%grass%daygrowth    = 0
        state%crop%grass%daygrowthpot = 0

! ---   actual rooting depth
        if (state%crop%common%swrd.eq.1) then
          state%crop%common%rd = afgen (state%crop%common%rdtb,22,rid)
          state%crop%common%rd = min(state%crop%common%rd,state%crop%common%rdm)
        elseif (state%crop%common%swrd.eq.2) then
          state%crop%common%rd = min(state%crop%common%rdi,state%crop%common%rdm)
        elseif (state%crop%common%swrd.eq.3) then
          state%crop%common%rdi = afgen (state%crop%common%rlwtb,22,state%crop%wofost%wrt)
          state%crop%common%rd = min(state%crop%common%rdi,state%crop%common%rdm)
        endif
        state%crop%common%rdpot = state%crop%common%rd
        
! ---   initial summation variables of the crop
        state%crop%wofost%tagp = state%crop%wofost%wlv+state%crop%wofost%wst
        state%crop%wofost%tagppot = state%crop%wofost%tagp
        state%crop%wofost%tagpt = 0.0d0
        state%crop%wofost%tagptpot = 0.0d0
        state%crop%common%cuptgraz = 0.0d0
        state%crop%common%cuptgrazpot = 0.0d0
        state%crop%common%tsum = 0.0d0
        
        state%crop%grass%cropstartpot     = rid
        state%crop%grass%cropstartact     = rid
        state%crop%grass%flhrvendpot      = .false.
        flearlyhrvendpot = .false.
        ! [SS-GR-CROP A5.1] mirror grass init-time state
        state%crop%common%cuptgraz    = state%crop%common%cuptgraz
        state%crop%common%cuptgrazpot = state%crop%common%cuptgrazpot
        
        if (state%crop%grass%swtsum.eq.0) then
          flGrassGrowth = .true.
        else
          flGrassGrowth = .false.  
        endif
        if (state%crop%grass%swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          call sumttd('initial',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif

! --- end skip above initialization if crop parameters are read from *.END file
      endif

      if (state%crop%swcf.ne.3) then
        state%crop%common%cf = afgen (state%crop%fixed%cftb,(2*magrs),rid)
        state%crop%common%ch = afgen (state%crop%fixed%chtb,(2*magrs),rid)
      else
        state%crop%common%cf        = afgen (state%crop%fixed%cftb,(2*magrs),state%crop%lai)
        state%crop%fixed%cfeic = afgen (state%crop%fixed%cfeictb,(2*magrs),state%crop%lai)
        state%crop%common%ch        = afgen(state%crop%fixed%chtb,(2*magrs),state%crop%lai)
      endif
      ! cfeic write retired (was: state%crop%fixed%cfeic = cfeic)

! --- initial storage on canopy
      if (state%crop%common%swinter.eq.3) then
        state%atmosphere%siccapact = siccaplai*state%crop%lai
      endif

! --- initialize matric flux potential (SS-CRP C-2.5: hroot/hleaf/mfluxtable
!     init moved to CropGrowth dispatcher which has access to state).
      if (state%crop%common%swdrought .eq. 2) then
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
        state%crop%grass%iharvest = 1
        do while (tc_t1900 .gt. state%crop%grass%dateharvest(state%crop%grass%iharvest))
          state%crop%grass%iharvest = state%crop%grass%iharvest + 1
        enddo
      endif      
      
! --- Find node for monitoring work-ability
      if (swlossmow .eq. 1) then
         
        ! Find node for monitoring work-ability during mowing
        nodmow = 1
        drz1       = -1.d0 * state%crop%grass%zmow - state%mesh%dz(nodmow)  ! [GR-BH C7]
        do while (drz1 .gt. 0.d0)
          nodmow   = nodmow + 1
          drz1 = drz1 - state%mesh%dz(nodmow)  ! [GR-BH C7]
        enddo
      
      endif   
      
      if (swlossgrz .eq. 1) then
        
        ! Find node and layer for monitoring work-ability at start of grazing
        nodgrz = 1
        drz1       = -1.d0 * state%crop%grass%zgrz - state%mesh%dz(nodgrz)  ! [GR-BH C7]
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
      if (state%crop%grass%flhrvendpot) then
        if (flearlyhrvendpot) then
          state%crop%grass%cropstartpot  = rid - 1.d0
        else
          state%crop%grass%cropstartpot  = rid
        endif
        state%crop%grass%pmowdm        = 0.d0
        state%crop%grass%pgrzdm        = 0.d0
        state%crop%wofost%plossdm       = 0.d0
      endif
      state%crop%grass%flhrvendpot      = .false.
      flearlyhrvendpot = .false.

! --- grass growth initiated by tsum from 1st day of calendar year
      state%crop%common%tsum = state%crop%common%tsum + max(0.0d0,at_tav)  ! [SS-GR-ATM B.5]
      if (.not. flGrassGrowth) then
        
        ! grass growth initiated by tsum
        if (state%crop%grass%swtsum.eq.1) then
          if (state%crop%common%tsum.ge.200.d0) then
            flGrassGrowth = .true.
          endif
        endif
        
        ! grass growth initiated by temperature, time and depth
        if (state%crop%grass%swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          if (dateGrassGrowth.eq.'undefined') call sumttd('dynamic',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif
      
        ! check if grass growth has started
        if (flGrassGrowth) then
          state%crop%grass%cropstartpot = rid
          state%crop%grass%cropstartact = rid
        endif

      endif
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.state%crop%grass%idregrpot) then

! ===   daily dry matter production ===

        gasspot = state%crop%wofost%pgasspot

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(state%crop%common%rmr*state%crop%wofost%wrtpot+state%crop%common%rml*state%crop%wofost%wlvpot+state%crop%common%rms*state%crop%wofost%wstpot)*afgen(state%crop%common%rfsetb,30,rid)
        teff = state%crop%common%q10**((at_tav-25.0d0)/10.0d0)  ! [SS-GR-ATM B.5]
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(state%crop%common%frtb,30,rid)
        fl = afgen(state%crop%common%fltb,30,rid)
        fs = afgen(state%crop%common%fstb,30,rid)
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
        cvf = 1.0d0/((fl/state%crop%common%cvl+fs/state%crop%common%cvs)*(1.0d0-fr)+fr/state%crop%common%cvr)
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
        if (state%crop%common%swrd.eq.3 .and. state%crop%wofost%wrtpot.gt.state%crop%common%wrtmax) then
          drrtpot = grrtpot
          drrtpot = max(drrtpot,state%crop%wofost%wrtpot*afgen (state%crop%common%rdrrtb,30,rid))
        else  
          drrtpot = state%crop%wofost%wrtpot*afgen (state%crop%common%rdrrtb,30,rid)
        endif  
        gwrtpot = grrtpot - drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/state%crop%kdif
        dslv2pot=state%crop%wofost%wlvpot*max(0.0d0,                                      &
     &                  min(0.03d0,0.03d0*(state%crop%common%laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = state%crop%common%ilvoldpot

        do while (restpot.gt.state%crop%common%lvpot(max(i1,1)).and.i1.ge.1)
          restpot = restpot-state%crop%common%lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (state%crop%common%lvagepot(max(i1,1)).gt.state%crop%common%span.and.restpot.gt.0.and.           &
     &                          i1.ge.1) then
          dalvpot = state%crop%common%lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.state%crop%common%lvagepot(max(i1,1)).gt.state%crop%common%span)
          dalvpot = dalvpot+state%crop%common%lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot   = dslvpot+dalvpot

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (state%crop%common%slatb,30,rid)
        if (state%crop%common%laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,at_tav-state%crop%common%tbase)  ! [SS-GR-ATM B.5]
          state%crop%common%glaiexpot = state%crop%common%laiexppot*state%crop%common%rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (state%crop%common%glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (state%crop%common%rdrstb,30,rid)*state%crop%wofost%wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        state%crop%grass%daygrowthpot = state%crop%grass%daygrowthpot + 1

!       Check trigger to start mowing event
        if (state%crop%grass%seqgrazmowpot(state%crop%grass%iseqgmpot) .eq. 2) then

          state%crop%grass%flHarvestpot = .false.
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (state%crop%wofost%tagppot .gt. dmharvest .or. (tc_daynr .gt. daylastharvest  &
     &          .and. state%crop%wofost%tagppot .gt. dmlastharvest)) then
                state%crop%grass%flHarvestpot = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(state%crop%grass%dmmowtb,20,rid)
              if (state%crop%wofost%tagppot .gt. dmharvest .or.                           &
     &                 (state%crop%grass%daygrowthpot .gt. maxdaymow .and. state%crop%grass%iseqgmpot .gt. 1)) then
                state%crop%grass%flHarvestpot = .true.
              endif
            endif

          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(tc_t1900 .gt. state%crop%grass%dateharvest(state%crop%grass%iharvest)) then
              state%crop%grass%flHarvestpot = .true.
            endif
          endif

!         In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
          if (state%crop%grass%flHarvestpot) then
            state%crop%grass%iseqgmpot = state%crop%grass%iseqgmpot + 1
            state%crop%common%slapot(1) = afgen (state%crop%common%slatb,30,rid)
            fl = afgen (state%crop%common%fltb,30,rid)
            fs = afgen (state%crop%common%fstb,30,rid)
            state%crop%wofost%wlvpot = state%crop%grass%mowrest / (1.d0 + (fs/fl))
            state%crop%wofost%wstpot = fs/fl*state%crop%wofost%wlvpot
            state%crop%wofost%dwlvpot = 0.0d0
            state%crop%wofost%dwstpot = 0.0d0
            state%crop%common%lvagepot(1) = 0.0d0
            state%crop%common%ilvoldpot = 1
            lasumpot = state%crop%wofost%wlvpot * state%crop%common%slapot(1)
            state%crop%common%laiexppot = lasumpot
            state%crop%common%lvpot(1) = state%crop%wofost%wlvpot
            
            gwstpot = 0.0d0
            gwrtpot = 0.0d0
            drlvpot = 0.0d0
            drstpot = 0.0d0
            drrtpot = 0.0d0
            
            state%crop%grass%daygrowthpot = 0
            
!           losses due to treading
            fralossmow = 0.d0
            if (swlossmow.eq.1) then
              fralossmow = afgen(state%crop%grass%lossmowtab,200,state%soilwater%h(nodmow))  ! [SS-SWC S-2.7]
            end if

!           harvest
            tagpspot = max(0.0d0,(state%crop%wofost%tagppot-(state%crop%wofost%wlvpot+state%crop%wofost%dwlvpot+state%crop%wofost%wstpot+state%crop%wofost%dwstpot)))
            state%crop%wofost%tagptpot = state%crop%wofost%tagptpot + tagpspot * (1.d0 - FraLossMow)

            state%crop%grass%cropendpot  = rid
            state%crop%grass%flhrvendpot = .true.
            state%crop%grass%pmowdm      = tagpspot * (1.d0 - FraLossMow)
            state%crop%wofost%plossdm     = tagpspot * FraLossMow
            
!           set regrowth delay
            idelaypot = int(afgen(state%crop%grass%DelayRegrowthTab,200,tagpspot))
            state%crop%grass%idregrpot = daycrop + idelaypot

          endif          
          
!       Check trigger to start grazing event          
        else if (state%crop%grass%seqgrazmowpot(state%crop%grass%iseqgmpot) .eq. 1 .or. state%crop%grass%seqgrazmowpot(state%crop%grass%iseqgmpot) .eq. 3) then

          state%crop%grass%flHarvestpot = .false.
          
          if (.not. state%crop%grass%flGrazingpot) then
              
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (state%crop%wofost%tagppot .gt. dmgrazing) then
                  state%crop%grass%flHarvestpot = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(state%crop%grass%dmgrztb,20,rid)
                if (state%crop%wofost%tagppot .gt. dmgrazing .or.                           &
     &                 (state%crop%grass%daygrowthpot .gt. maxdaygrz .and. state%crop%grass%iseqgmpot .gt. 1)) then
                  state%crop%grass%flHarvestpot = .true.
                endif
              endif
              
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(tc_t1900 .gt. state%crop%grass%dateharvest(state%crop%grass%iharvest)) then
                state%crop%grass%flHarvestpot = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (state%crop%grass%flHarvestpot .or. state%crop%grass%flGrazingpot) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgrazpot = state%crop%grass%lsda(state%crop%grass%iseqgmpot) *                                 &
     &                         afgen(state%crop%grass%uptgrazingtab,200,state%crop%grass%lsda(state%crop%grass%iseqgmpot))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgrazpot = state%crop%grass%lsda(state%crop%grass%iseqgmpot) *                                &
     &                         afgen(state%crop%grass%lossgrazingtab,200,state%crop%grass%lsda(state%crop%grass%iseqgmpot))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(state%crop%grass%lossgrztab,200,state%soilwater%h(nodgrz))  ! [SS-SWC S-2.7]
            end if
            lossgrazpot = lossgrazpot + state%crop%wofost%tagppot * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. state%crop%grass%flGrazingpot) then
              state%crop%grass%daygrowthpot   = 0
              state%crop%grass%idaysgrazpot   = 0
              flDewoolingpot = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((state%crop%wofost%tagppot - uptgrazpot - lossgrazpot) .gt. tagprest) then
              
              state%crop%grass%flGrazingpot = .true.
              state%crop%common%cuptgrazpot  = state%crop%common%cuptgrazpot + uptgrazpot
          
!             distribute grazing over stems and leaves (living and dead parts)
              state%crop%wofost%wstpot  = state%crop%wofost%wstpot  - (uptgrazpot+lossgrazpot) * state%crop%wofost%wstpot  / state%crop%wofost%tagppot
              state%crop%wofost%dwstpot = state%crop%wofost%dwstpot - (uptgrazpot+lossgrazpot) * state%crop%wofost%dwstpot / state%crop%wofost%tagppot
              state%crop%wofost%dwlvpot = state%crop%wofost%dwlvpot - (uptgrazpot+lossgrazpot) * state%crop%wofost%dwlvpot / state%crop%wofost%tagppot
              grazlivinglvpot =   (uptgrazpot+lossgrazpot) * state%crop%wofost%wlvpot  / state%crop%wofost%tagppot
          
!             reduce leave weights
              i1 = state%crop%common%ilvoldpot
              do while (grazlivinglvpot .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglvpot .ge. state%crop%common%lvpot(i1)) then
                  grazlivinglvpot = grazlivinglvpot - state%crop%common%lvpot(i1)
                  state%crop%common%lvpot(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  state%crop%common%lvpot(i1) = state%crop%common%lvpot(i1) - grazlivinglvpot
                  grazlivinglvpot = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              state%crop%grass%cropendpot = rid
              state%crop%grass%pgrzdm     = state%crop%grass%pgrzdm + uptgrazpot
              state%crop%wofost%plossdm    = state%crop%wofost%tagppot * fralossgrz
              
!             Check number of days with grazing
              daysgrazpot  = int(afgen(state%crop%grass%daysgrazingtab,200,state%crop%grass%lsda(state%crop%grass%iseqgmpot)))
              state%crop%grass%idaysgrazpot = state%crop%grass%idaysgrazpot + 1
              if(state%crop%grass%idaysgrazpot .eq. daysgrazpot) then
                state%crop%grass%flGrazingpot = .false.
                state%crop%grass%flhrvendpot  = .true.
                if (state%crop%grass%seqgrazmowpot(state%crop%grass%iseqgmpot) .eq. 3) then
                  flDewoolingpot  = .true.
                endif
                state%crop%grass%daygrowthpot = 0
                state%crop%grass%iseqgmpot = state%crop%grass%iseqgmpot + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (state%crop%grass%flGrazingpot .or. swharvest .eq. 2) then
              state%crop%grass%flGrazingpot     = .false.
              state%crop%grass%flhrvendpot      = .true.
              flearlyhrvendpot = .true.
              if (state%crop%grass%seqgrazmowpot(state%crop%grass%iseqgmpot) .eq. 3 .and. state%crop%wofost%tagppot .gt. state%crop%grass%dewrest) then
                flDewoolingpot   = .true.
                flearlyhrvendpot = .false.
              endif
              state%crop%grass%daygrowthpot = 0
              state%crop%grass%iseqgmpot = state%crop%grass%iseqgmpot + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            state%crop%grass%idregrpot = daycrop

!           Dewooling after grazing event            
            if (flDewoolingpot) then

              state%crop%common%slapot(1) = afgen (state%crop%common%slatb,30,rid)
              fl = afgen (state%crop%common%fltb,30,rid)
              fs = afgen (state%crop%common%fstb,30,rid)
              state%crop%wofost%wlvpot = state%crop%grass%dewrest / (1.d0 + (fs/fl))
              state%crop%wofost%wstpot = fs/fl*state%crop%wofost%wlvpot
              state%crop%wofost%dwlvpot = 0.0d0
              state%crop%wofost%dwstpot = 0.0d0
              state%crop%common%lvagepot(1) = 0.0d0
              state%crop%common%ilvoldpot = 1
              lasumpot = state%crop%wofost%wlvpot * state%crop%common%slapot(1)
              state%crop%common%laiexppot = lasumpot
              state%crop%common%lvpot(1) = state%crop%wofost%wlvpot
              
              gwstpot = 0.0d0
              gwrtpot = 0.0d0
              drlvpot = 0.0d0
              drstpot = 0.0d0
              drrtpot = 0.0d0
              
!             Assumption: one day delay in regrowth after grazing
              state%crop%grass%idregrpot = daycrop + 1
              
            endif
            
          endif
          
        endif
        
        if (daycrop .ge. state%crop%grass%idregrpot) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(at_tav-state%crop%common%tbase)/(35.0d0-state%crop%common%tbase))  ! [SS-GR-ATM B.5]

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = state%crop%common%ilvoldpot
           do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.state%crop%common%lvpot(i1)) then
              dslvtpot = dslvtpot-state%crop%common%lvpot(i1)
              state%crop%common%lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              state%crop%common%lvpot(i1) = state%crop%common%lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (state%crop%common%lvagepot(max(i1,1)) .gt. state%crop%common%span .and. i1 .ge. 1)
              state%crop%common%lvpot(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          state%crop%common%ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = state%crop%common%ilvoldpot,1,-1
            state%crop%common%lvpot(i1+1) = state%crop%common%lvpot(i1)
            state%crop%common%slapot(i1+1) = state%crop%common%slapot(i1)
            state%crop%common%lvagepot(i1+1) = state%crop%common%lvagepot(i1)+fysdel*delt
          enddo
          state%crop%common%ilvoldpot = state%crop%common%ilvoldpot+1

! ---     new leaves in class 1
          state%crop%common%lvpot(1) = grlvpot*delt
          state%crop%common%slapot(1) = slatpot
          state%crop%common%lvagepot(1) = 0.d0

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          state%crop%wofost%wlvpot = 0.d0
          do i1 = 1,state%crop%common%ilvoldpot
            lasumpot = lasumpot+state%crop%common%lvpot(i1)*state%crop%common%slapot(i1)
            state%crop%wofost%wlvpot = state%crop%wofost%wlvpot+state%crop%common%lvpot(i1)
          enddo

          state%crop%common%laiexppot = state%crop%common%laiexppot+state%crop%common%glaiexpot*delt

        endif

! ---   dry weight of living plant organs
        state%crop%wofost%wrtpot = state%crop%wofost%wrtpot+gwrtpot*delt
        state%crop%wofost%wstpot = state%crop%wofost%wstpot+gwstpot*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        state%crop%wofost%dwrtpot = state%crop%wofost%dwrtpot+drrtpot*delt
        state%crop%wofost%dwlvpot = state%crop%wofost%dwlvpot+drlvpot*delt
        state%crop%wofost%dwstpot = state%crop%wofost%dwstpot+drstpot*delt

! ---   dry weight of dead and living plant organs
        twlvpot = state%crop%wofost%wlvpot+state%crop%wofost%dwlvpot
        twstpot = state%crop%wofost%wstpot+state%crop%wofost%dwstpot
        state%crop%wofost%tagppot = twlvpot+twstpot

! ---   leaf area index
        state%crop%common%laipot = lasumpot+state%crop%common%ssa*state%crop%wofost%wstpot
!       prevent immediate lai reduction at emergence
!       KRO-BOO-20160403: suppressed because deviates from Wofost
!       laipot = max(laipot, laiem)

        ! root extension
        if (state%crop%common%swrd.eq.1) then
          state%crop%common%rdpot = afgen (state%crop%common%rdtb,22,rid)
          state%crop%common%rdpot = min(state%crop%common%rdpot,state%crop%common%rdm)
        elseif (state%crop%common%swrd.eq.2) then
          rrpot = min (state%crop%common%rdm-state%crop%common%rdpot,state%crop%common%rri)
          if (fr.le.0.0d0 .or. state%crop%wofost%pgasspot.lt.1.0d0) rrpot = 0.0d0
          state%crop%common%rdpot = state%crop%common%rdpot + rrpot
        elseif (state%crop%common%swrd.eq.3) then
          state%crop%common%rdpot = afgen (state%crop%common%rlwtb,22,state%crop%wofost%wrtpot)
          state%crop%common%rdpot = min(state%crop%common%rdpot,state%crop%common%rdm)
        endif

      endif

      ! [SS-GR-CROP A5.1] mirror grass case(2) potential state
      state%crop%common%cuptgrazpot   = state%crop%common%cuptgrazpot

      return

      case (3)

! === calculate actual rate and state variables ======================================

! --- check end of harvest
      if (state%crop%grass%flhrvendact) then
        if (flearlyhrvendact) then
          state%crop%grass%cropstartact  = rid - 1.d0
        else
          state%crop%grass%cropstartact  = rid
        endif
        state%crop%grass%mowdm        = 0.d0
        state%crop%grass%grzdm        = 0.d0
        state%crop%wofost%lossdm       = 0.d0
      endif
      state%crop%grass%flhrvendact      = .false.
      flearlyhrvendact = .false.

! --- rates of change of the crop variables ---------------------------------------------
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop .ge. state%crop%grass%idregr) then

! ===   daily dry matter production ===

! ---   water stress reduction of pgass to gass
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if(dabs(state%atmosphere%ptra).lt.nihil) then
          reltr = 1.0d0
        else
          reltr = max(0.0d0,min(1.0d0,state%soilwater%tra/state%atmosphere%ptra))  ! [SS-SWC S-2.7]
        endif
        gass = state%crop%wofost%pgass * reltr

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (state%crop%common%rmr*state%crop%wofost%wrt+state%crop%common%rml*state%crop%wofost%wlv+state%crop%common%rms*state%crop%wofost%wst)*afgen(state%crop%common%rfsetb,30,rid)
        teff = state%crop%common%q10**((at_tav-25.0d0)/10.0d0)  ! [SS-GR-ATM B.5]
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   partitioning factors (relevant for restart)
        fr = afgen(state%crop%common%frtb,30,rid)
        fl = afgen(state%crop%common%fltb,30,rid)
        fs = afgen(state%crop%common%fstb,30,rid)

! ---   dry matter increase
        cvf = 1.0d0/((fl/state%crop%common%cvl+fs/state%crop%common%cvs)*(1.0d0-fr)+fr/state%crop%common%cvr)
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
        if (state%crop%common%swrd.eq.3 .and. state%soilwater%flWrtNonox) grrt = 0.d0   ! [SS-GR-CROPWS A4] present(state) guard removed
        if (state%crop%common%swrd.eq.3 .and. state%crop%wofost%wrt.gt.state%crop%common%wrtmax) then
          drrt = grrt
          drrt = max(drrt,state%crop%wofost%wrt*afgen (state%crop%common%rdrrtb,30,rid))
        else  
          drrt = state%crop%wofost%wrt*afgen (state%crop%common%rdrrtb,30,rid)
        endif  
        gwrt = grrt-drrt

! ---   growth rate leaves

! ---   weight of new leaves
        admi = (1.0d0-fr)*dmi        
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        dslv1 = state%crop%wofost%wlv*(1.0d0-reltr)*state%crop%common%perdl
        laicr = 3.2d0/state%crop%kdif
        dslv2 = state%crop%wofost%wlv*max(0.0d0,min(0.03d0,0.03d0*(state%crop%lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = state%crop%common%ilvold

        do while (rest.gt.state%crop%common%lv(max(i1,1)).and.i1.ge.1)
          rest = rest-state%crop%common%lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (state%crop%common%lvage(max(i1,1)).gt.state%crop%common%span.and.rest.gt.0.and.i1.ge.1) then
          dalv = state%crop%common%lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.state%crop%common%lvage(max(i1,1)).gt.state%crop%common%span)
          dalv = dalv+state%crop%common%lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv   = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (state%crop%common%slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (state%crop%common%laiexp.lt.6.0d0) then
          dteff = max (0.0d0,at_tav-state%crop%common%tbase)  ! [SS-GR-ATM B.5]
          state%crop%common%glaiex = state%crop%common%laiexp*state%crop%common%rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (state%crop%common%glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        drst1 = state%crop%wofost%wst*(1.0d0-reltr)*state%crop%common%perdl
! ---   death of stems due to ageing
        drst2 = afgen (state%crop%common%rdrstb,30,rid)*state%crop%wofost%wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        state%crop%grass%daygrowth = state%crop%grass%daygrowth + 1

!       Check trigger to start mowing event
        if (state%crop%grass%seqgrazmow(state%crop%grass%iseqgm) .eq. 2) then
          
          state%crop%grass%flHarvest = .false.   
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (state%crop%wofost%tagp .gt. dmharvest .or. (tc_daynr .gt. daylastharvest  &
     &          .and. state%crop%wofost%tagp .gt. dmlastharvest)) then
                state%crop%grass%flHarvest = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(state%crop%grass%dmmowtb,20,rid)
              if (state%crop%wofost%tagp .gt. dmharvest .or.                           &
     &                 (state%crop%grass%daygrowth .gt. maxdaymow .and. state%crop%grass%iseqgm .gt. 1)) then
                state%crop%grass%flHarvest = .true.
              endif
            endif
          
          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(tc_t1900 .gt. state%crop%grass%dateharvest(state%crop%grass%iharvest)) then
              state%crop%grass%iharvest = state%crop%grass%iharvest + 1
              state%crop%grass%flHarvest = .true.
            endif
          endif
          
!       In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
        if (state%crop%grass%flHarvest) then
          state%crop%grass%iseqgm = state%crop%grass%iseqgm + 1
          state%crop%common%sla(1) = afgen (state%crop%common%slatb,30,rid)
          fl = afgen (state%crop%common%fltb,30,rid)
          fs = afgen (state%crop%common%fstb,30,rid)
          state%crop%wofost%wlv = state%crop%grass%mowrest / (1.d0 + (fs/fl))
          state%crop%wofost%wst = fs/fl*state%crop%wofost%wlv
          state%crop%wofost%dwlv = 0.0d0
          state%crop%wofost%dwst = 0.0d0
          state%crop%common%lvage(1) = 0.0d0
          state%crop%common%ilvold = 1
          lasum = state%crop%wofost%wlv * state%crop%common%sla(1)
          state%crop%common%laiexp = lasum
          state%crop%common%lv(1) = state%crop%wofost%wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          state%crop%grass%daygrowth = 0

!         losses due to treading
          FraLossMow = 0.d0
          if (swlossmow.eq.1) then
            FraLossMow = afgen(state%crop%grass%lossmowtab,200,state%soilwater%h(nodmow))  ! [SS-SWC S-2.7]
          end if
          
!         harvest
          tagps = max (0.0d0,(state%crop%wofost%tagp-(state%crop%wofost%wlv+state%crop%wofost%dwlv+state%crop%wofost%wst+state%crop%wofost%dwst)))
          state%crop%wofost%tagpt = state%crop%wofost%tagpt + tagps * (1.d0 - fralossmow)

          state%crop%grass%cropendact  = rid
          state%crop%grass%flhrvendact = .true.
          state%crop%grass%mowdm   = tagps * (1.d0 - FraLossMow)
          state%crop%wofost%lossdm  = tagps * FraLossMow
          
! ---     set regrowth delay
          idelay = int(afgen(state%crop%grass%DelayRegrowthTab,200,tagps))
          state%crop%grass%idregr = daycrop + idelay

        endif
          
!       Check trigger to start grazing event          
        else if (state%crop%grass%seqgrazmow(state%crop%grass%iseqgm) .eq. 1 .or. state%crop%grass%seqgrazmow(state%crop%grass%iseqgm) .eq. 3) then

          state%crop%grass%flHarvest = .false.
            
          if (.not. state%crop%grass%flGrazing) then
            
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (state%crop%wofost%tagp .gt. dmgrazing) then
                  state%crop%grass%flHarvest = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(state%crop%grass%dmgrztb,20,rid)
                if (state%crop%wofost%tagp .gt. dmgrazing .or.                           &
     &            (state%crop%grass%daygrowth .gt. maxdaygrz .and. state%crop%grass%iseqgm .gt. 1)) then
                  state%crop%grass%flHarvest = .true.
                endif
              endif
            
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(tc_t1900 .gt. state%crop%grass%dateharvest(state%crop%grass%iharvest)) then
                state%crop%grass%iharvest = state%crop%grass%iharvest + 1
                state%crop%grass%flHarvest = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (state%crop%grass%flHarvest .or. state%crop%grass%flGrazing) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = state%crop%grass%lsda(state%crop%grass%iseqgm)*afgen(state%crop%grass%uptgrazingtab,200,state%crop%grass%lsda(state%crop%grass%iseqgm))            

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = state%crop%grass%lsda(state%crop%grass%iseqgm) *                                &
     &                         afgen(state%crop%grass%lossgrazingtab,200,state%crop%grass%lsda(state%crop%grass%iseqgm))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(state%crop%grass%lossgrztab,200,state%soilwater%h(nodgrz))  ! [SS-SWC S-2.7]
            end if
            lossgraz = lossgraz + state%crop%wofost%tagp * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. state%crop%grass%flGrazing) then
              state%crop%grass%daygrowth   = 0
              state%crop%grass%idaysgraz   = 0
              flDewooling = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((state%crop%wofost%tagp - uptgraz - lossgraz) .gt. tagprest) then
              
              state%crop%grass%flGrazing = .true.
              state%crop%common%cuptgraz  = state%crop%common%cuptgraz + uptgraz

!             distribute grazing over stems and leaves (living and dead parts)
              state%crop%wofost%wst  = state%crop%wofost%wst  -  (uptgraz+lossgraz) * state%crop%wofost%wst  / state%crop%wofost%tagp
              state%crop%wofost%dwst = state%crop%wofost%dwst -  (uptgraz+lossgraz) * state%crop%wofost%dwst / state%crop%wofost%tagp
              state%crop%wofost%dwlv = state%crop%wofost%dwlv -  (uptgraz+lossgraz) * state%crop%wofost%dwlv / state%crop%wofost%tagp
              grazlivinglv = (uptgraz+lossgraz) * state%crop%wofost%wlv  / state%crop%wofost%tagp
          
!             reduce leave weights
              i1 = state%crop%common%ilvold
              do while (grazlivinglv .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglv .ge. state%crop%common%lv(i1)) then
                  grazlivinglv = grazlivinglv - state%crop%common%lv(i1)
                  state%crop%common%lv(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  state%crop%common%lv(i1) = state%crop%common%lv(i1) - grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              state%crop%grass%cropendact = rid
              state%crop%grass%grzdm      = state%crop%grass%grzdm + uptgraz
              state%crop%wofost%lossdm     = state%crop%wofost%tagp * fralossgrz
              
!             Check number of days with grazing
              daysgraz  = int(afgen(state%crop%grass%daysgrazingtab,200,state%crop%grass%lsda(state%crop%grass%iseqgm)))
              state%crop%grass%idaysgraz = state%crop%grass%idaysgraz + 1
              if(state%crop%grass%idaysgraz .eq. daysgraz) then
                state%crop%grass%flGrazing   = .false.
                state%crop%grass%flhrvendact = .true.
                if (state%crop%grass%seqgrazmow(state%crop%grass%iseqgm) .eq. 3) then
                  flDewooling  = .true.
                endif
                state%crop%grass%daygrowth = 0
                state%crop%grass%iseqgm = state%crop%grass%iseqgm + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (state%crop%grass%flGrazing .or. swharvest .eq. 2) then
              state%crop%grass%flGrazing        = .false.
              state%crop%grass%flhrvendact      = .true.
              flearlyhrvendact = .true.
              if (state%crop%grass%seqgrazmow(state%crop%grass%iseqgm) .eq. 3 .and. state%crop%wofost%tagp .gt. state%crop%grass%dewrest) then
                flDewooling      = .true.
                flearlyhrvendact = .false.
              endif
              state%crop%grass%daygrowth = 0
              state%crop%grass%iseqgm = state%crop%grass%iseqgm + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            state%crop%grass%idregr = daycrop

!           Dewooling after grazing event            
            if (flDewooling) then

              state%crop%common%sla(1) = afgen (state%crop%common%slatb,30,rid)
              fl = afgen (state%crop%common%fltb,30,rid)
              fs = afgen (state%crop%common%fstb,30,rid)
              state%crop%wofost%wlv = state%crop%grass%dewrest / (1.d0 + (fs/fl))
              state%crop%wofost%wst = fs/fl*state%crop%wofost%wlv
              state%crop%wofost%dwlv = 0.0d0
              state%crop%wofost%dwst = 0.0d0
              state%crop%common%lvage(1) = 0.0d0
              state%crop%common%ilvold = 1
              lasum = state%crop%wofost%wlv * state%crop%common%sla(1)
              state%crop%common%laiexp = lasum
              state%crop%common%lv(1) = state%crop%wofost%wlv
    
              gwst = 0.0d0
              gwrt = 0.0d0
              drlv = 0.0d0
              drst = 0.0d0
              drrt = 0.0d0
    
!             Assumption: one day delay in regrowth after grazing
              state%crop%grass%idregr = daycrop + 1

            endif
            
          endif
          
        endif

        if (daycrop .ge. state%crop%grass%idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(at_tav-state%crop%common%tbase)/(35.0d0-state%crop%common%tbase))  ! [SS-GR-ATM B.5]

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = state%crop%common%ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.state%crop%common%lv(i1)) then
              dslvt = dslvt-state%crop%common%lv(i1)
              state%crop%common%lv(i1) = 0.0d0
              i1 = i1-1
            else
              state%crop%common%lv(i1) = state%crop%common%lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (state%crop%common%lvage(max(i1,1)).gt.state%crop%common%span.and.i1.ge.1)
              state%crop%common%lv(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          state%crop%common%ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = state%crop%common%ilvold,1,-1
            state%crop%common%lv(i1+1) = state%crop%common%lv(i1)
            state%crop%common%sla(i1+1) = state%crop%common%sla(i1)
            state%crop%common%lvage(i1+1) = state%crop%common%lvage(i1)+fysdel*delt
          enddo
          state%crop%common%ilvold = state%crop%common%ilvold+1

! ---     new leaves in class 1
          state%crop%common%lv(1) = grlv*delt
          state%crop%common%sla(1) = slat
          state%crop%common%lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          state%crop%wofost%wlv = 0.d0
          do i1 = 1,state%crop%common%ilvold
            lasum = lasum+state%crop%common%lv(i1)*state%crop%common%sla(i1)
            state%crop%wofost%wlv = state%crop%wofost%wlv+state%crop%common%lv(i1)
          enddo

          state%crop%common%laiexp = state%crop%common%laiexp+state%crop%common%glaiex*delt

        endif

! ---   dry weight of living plant organs
        state%crop%wofost%wrt = state%crop%wofost%wrt+gwrt*delt
        state%crop%wofost%wst = state%crop%wofost%wst+gwst*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        state%crop%wofost%dwrt = state%crop%wofost%dwrt+drrt*delt
        state%crop%wofost%dwlv = state%crop%wofost%dwlv+drlv*delt
        state%crop%wofost%dwst = state%crop%wofost%dwst+drst*delt

! ---   dry weight of dead and living plant organs
        twlv = state%crop%wofost%wlv+state%crop%wofost%dwlv
        twst = state%crop%wofost%wst+state%crop%wofost%dwst
        state%crop%wofost%tagp = twlv+twst

! ---   leaf area index
        state%crop%lai = lasum+state%crop%common%ssa*state%crop%wofost%wst
        state%crop%common%laimax = max (state%crop%lai,state%crop%common%laimax)

! ---   update normalized cumulative root density based on root extraction or stress (cumdens)
        if (state%crop%common%swrdc .eq. 1) call update_rootdistribution(state)
        
        ! root extension
        if (state%crop%common%swrd.eq.1) then
          state%crop%common%rd = afgen (state%crop%common%rdtb,22,rid)
          state%crop%common%rd = min(state%crop%common%rd,state%crop%common%rdm)
        elseif (state%crop%common%swrd.eq.2) then
          rr = min (state%crop%common%rdm-state%crop%common%rd,state%crop%common%rri)
          if (fr.le.0.0d0 .or. state%crop%wofost%pgass.lt.1.0d0 .or.                    &
     &        state%soilwater%flWrtNonox) rr = 0.0d0   ! [SS-GR-CROPWS A4] present(state) guard removed
          if (state%crop%common%swdmi2rd.eq.1 .and. state%crop%wofost%pgass.ge.1.0d0)              rr = rr * gass/state%crop%wofost%pgass
          state%crop%common%rd = state%crop%common%rd + rr
        elseif (state%crop%common%swrd.eq.3) then
          state%crop%common%rd = afgen (state%crop%common%rlwtb,22,state%crop%wofost%wrt)
          state%crop%common%rd = min(state%crop%common%rd,state%crop%common%rdm)
        endif

! ---   set crop height and cropfactor
        if (state%crop%swcf.ne.3) then
          state%crop%common%cf = afgen (state%crop%fixed%cftb,(2*magrs),rid)
          state%crop%common%ch = afgen (state%crop%fixed%chtb,(2*magrs),rid)
        else
          state%crop%common%cf = afgen (state%crop%fixed%cftb,(2*magrs),state%crop%lai)
          state%crop%fixed%cfeic = afgen (state%crop%fixed%cfeictb,(2*magrs),state%crop%lai)
          state%crop%common%ch = afgen(state%crop%fixed%chtb,(2*magrs),state%crop%lai)
        endif
        ! cfeic write retired (was: state%crop%fixed%cfeic = cfeic)

! ---   update canopy storage capacity
        if (state%crop%common%swinter.eq.3) then
          state%atmosphere%siccapact = siccaplai*state%crop%lai
        endif

      endif

      ! [SS-GR-CROP A5.1] mirror grass case(3) actual state
      state%crop%wofost%tagp        = state%crop%wofost%tagp
      state%crop%common%cuptgraz    = state%crop%common%cuptgraz

      return

      case default
         call fatalerr_collected ('Grass', 'Illegal value for TASK')
      end select

      end associate  ! tc_t1900, tc_daynr => state%timecontrol [TC-10]
      return
      end subroutine grass

      end module cropgrass_runtime_mod
