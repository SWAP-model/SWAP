! cropgrass_runtime.f90 — type-3 (grass / WOFOST grass) crop runtime dispatcher.
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
!
! [GR-CROP 2026-05-25] crop-sweep:
!   - magrs/macp sourced from swap_array_dimensions.
!   - rdmax read via state%cfg%crop%rdmax (Class B direct read).
!   - siccaplai*lai branches: swinter=3 stub-errored on TOML; replaced with 0.0d0.
!   - wiltpoint → state%crop%common%hlim4 (Feddes wilting-point head).
!   - twilt/flhydrlift migrated to state%crop%common (drought-stress workspace).
!   - reltr migrated to state%crop%common%reltr (active rotation only — written
!     by the cropX_runtime task=3 selected by the dispatcher; read same-file).
!   - Remaining legacy globals (rid, daycrop, wrtmin, gwrt): cross-file readers
!     in oxygenstress / cropgrowth_helpers / cropwofost_runtime.
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: magrs, macp
      use variables, only: rid, daycrop,        &  ! workspace/scratch — cross-file readers (oxygenstress, cropgrowth)
                           wrtmin, gwrt           ! cross-file with cropgrowth_helpers/cropwofost_runtime
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      ! GR-CROPWS Phase 0: sumttd, update_rootdistribution extracted to cropgrowth_helpers_mod
      use cropgrowth_helpers_mod, only: sumttd, update_rootdistribution

      implicit none

      type(swap_state_t), intent(inout) :: state 

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
      ! [GR-CROP 2026-05-25] sub-record associate (crop-sweep convention).
      associate( crop     => state%crop,            &
                 soil     => state%soilwater,       &
                 mesh     => state%mesh,            &
                 atmo     => state%atmosphere,      &
                 time     => state%timecontrol,     &
                 cfg_crop => state%cfg%crop         )

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
               if (crop%common%icrop >= 1 .and. crop%common%icrop <= size(crop_config_global%rotation_loaded)) then
                  if (crop_config_global%rotation_loaded(crop%common%icrop)) then
                     ! Defense-in-depth: only dispatch to cache when the schema
                     ! is fully authored (case 4 + case 2 have amaxtb; the
                     ! hupselbrook skeleton does not — Phase 4 will fill it).
                     if (allocated(crop_config_global%rotation_grass(crop%common%icrop)%amaxtb)) &
                        use_cache = .true.
                  end if
               end if
            end if
         end if
         if (use_cache) then
            associate(cfg => crop_config_global%rotation_grass(crop%common%icrop))
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
               call cropgrass_init_from_config(cfg, crop%common%icrop, &
                  time%tend, time%tstart, state)
            end associate
         else
            ! ADR 0016 cache-miss: typed config required for type=3 rotations.
            ! No silent legacy fallback — the user must author cropgrass.crp.toml.
            call fatalerr_collected('cropgrowth/Grass', &
               'cropgrass rotation has no loaded .crp.toml — author the file or use the legacy executable.')
         end if
      end block

! --- sequence of harvest by mowing, dewooling and grazing
      crop%grass%seqgrazmowpot = crop%grass%seqgrazmow

! --- development stage (not used by Grassland, instead Daynrs are used)
      crop%common%dvs = -99.99d0

! --- maximum rooting depth
      if (crop%common%swrd.eq.1) then
        crop%common%rdm = cfg_crop%rdmax
      elseif (crop%common%swrd.eq.2) then
        crop%common%rdm = min(cfg_crop%rdmax,crop%common%rdc)
      elseif (crop%common%swrd.eq.3) then
        crop%common%rdc = afgen (crop%common%rlwtb,22,crop%common%wrtmax)
        crop%common%rdm = min(cfg_crop%rdmax,crop%common%rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (time%t1900 - time%tstart .gt. tiny .or. soil%swinco .ne. 3 .or.     &
     &   dabs(time%t1900 - crop%common%cropstart) .lt. tiny) then 

        crop%grass%iseqgm = 1
        crop%grass%iseqgmpot = crop%grass%iseqgm

! ---   initial values of crop parameters
        rid = dble(daycrop)
        fr = afgen (crop%common%frtb,30,rid)
        fl = afgen (crop%common%fltb,30,rid)
        fs = afgen (crop%common%fstb,30,rid)
        crop%common%sla(1) = afgen (crop%common%slatb,30,rid)
        crop%common%lvage(1) = 0.d0
        crop%common%ilvold = 1
        crop%grass%idregr = 0
        crop%common%slapot(1) = afgen (crop%common%slatb,30,rid)
        crop%common%lvagepot(1) = 0.d0
        crop%common%ilvoldpot = 1
        crop%grass%idregrpot = 0

! ---   initial state variables of the crop
        crop%wofost%wrt = fr*crop%common%tdwi
        wrtmin = crop%wofost%wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        crop%wofost%wrtpot = crop%wofost%wrt
        crop%wofost%wst = fs*(1.0d0-fr)*crop%common%tdwi
        crop%wofost%wstpot = crop%wofost%wst
        crop%wofost%wlv = crop%common%laiem/crop%common%sla(1)
        crop%wofost%wlvpot = crop%wofost%wlv
        
!     KRO-BOO-20160403: intro because comparison with Wofost
        crop%common%laiem = crop%wofost%wlv*crop%common%sla(1)  ! is not input !
        crop%common%lv(1) = crop%wofost%wlv
        crop%common%lvpot(1) = crop%common%lv(1)
        lasum = crop%common%laiem
        lasumpot = lasum     
        crop%common%glaiex = 0.0d0
        crop%common%glaiexpot = 0.0d0
        crop%common%laiexp = crop%common%laiem
        crop%common%laiexppot = crop%common%laiem
        crop%common%laimax = crop%common%laiem
        crop%lai = lasum+crop%common%ssa*crop%wofost%wst
        crop%common%laipot = crop%lai
        crop%wofost%dwrt = 0.d0
        crop%wofost%dwrtpot = crop%wofost%dwrt
        crop%wofost%dwlv = 0.d0
        crop%wofost%dwlvpot = crop%wofost%dwlv
        crop%wofost%dwst = 0.d0
        crop%wofost%dwstpot = crop%wofost%dwst

        crop%grass%daygrowth    = 0
        crop%grass%daygrowthpot = 0

! ---   actual rooting depth
        if (crop%common%swrd.eq.1) then
          crop%common%rd = afgen (crop%common%rdtb,22,rid)
          crop%common%rd = min(crop%common%rd,crop%common%rdm)
        elseif (crop%common%swrd.eq.2) then
          crop%common%rd = min(crop%common%rdi,crop%common%rdm)
        elseif (crop%common%swrd.eq.3) then
          crop%common%rdi = afgen (crop%common%rlwtb,22,crop%wofost%wrt)
          crop%common%rd = min(crop%common%rdi,crop%common%rdm)
        endif
        crop%common%rdpot = crop%common%rd
        
! ---   initial summation variables of the crop
        crop%wofost%tagp = crop%wofost%wlv+crop%wofost%wst
        crop%wofost%tagppot = crop%wofost%tagp
        crop%wofost%tagpt = 0.0d0
        crop%wofost%tagptpot = 0.0d0
        crop%common%cuptgraz = 0.0d0
        crop%common%cuptgrazpot = 0.0d0
        crop%common%tsum = 0.0d0
        
        crop%grass%cropstartpot     = rid
        crop%grass%cropstartact     = rid
        crop%grass%flhrvendpot      = .false.
        flearlyhrvendpot = .false.
        ! mirror grass init-time state
        crop%common%cuptgraz    = crop%common%cuptgraz
        crop%common%cuptgrazpot = crop%common%cuptgrazpot
        
        if (crop%grass%swtsum.eq.0) then
          flGrassGrowth = .true.
        else
          flGrassGrowth = .false.  
        endif
        if (crop%grass%swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          call sumttd('initial',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif

! --- end skip above initialization if crop parameters are read from *.END file
      endif

      if (crop%swcf.ne.3) then
        crop%common%cf = afgen (crop%fixed%cftb,(2*magrs),rid)
        crop%common%ch = afgen (crop%fixed%chtb,(2*magrs),rid)
      else
        crop%common%cf        = afgen (crop%fixed%cftb,(2*magrs),crop%lai)
        crop%fixed%cfeic = afgen (crop%fixed%cfeictb,(2*magrs),crop%lai)
        crop%common%ch        = afgen(crop%fixed%chtb,(2*magrs),crop%lai)
      endif
      ! cfeic write retired (was: crop%fixed%cfeic = cfeic)

! --- initial storage on canopy
      if (crop%common%swinter.eq.3) then
        atmo%siccapact = 0.0d0   ! [GR-CROP 2026-05-25] swinter=3 stub-errored on TOML; siccaplai always 0
      endif

! --- initialize matric flux potential (SS-CRP C-2.5: hroot/hleaf/mfluxtable
!     init moved to CropGrowth dispatcher which has access to state).
      if (crop%common%swdrought .eq. 2) then
        if (swhydrlift .eq. 1) then
          crop%common%flhydrlift = .true.
        else
          crop%common%flhydrlift = .false.
        endif
        if (.not. allocated(crop%common%twilt)) allocate(crop%common%twilt(mesh%numnod))
        do i = 1,mesh%numnod
         crop%common%twilt(i) = watcon(crop%common%hlim4, &
                            soil%vg_params(i), &
                            soil%iHWCKmodel(soil%layer(i)), &
                            i, soil)
        enddo
      endif

! --- harvest
!     initialise 
      if (swharvest.eq.2) then
        crop%grass%iharvest = 1
        do while (time%t1900 .gt. crop%grass%dateharvest(crop%grass%iharvest))
          crop%grass%iharvest = crop%grass%iharvest + 1
        enddo
      endif      
      
! --- Find node for monitoring work-ability
      if (swlossmow .eq. 1) then
         
        ! Find node for monitoring work-ability during mowing
        nodmow = 1
        drz1       = -1.d0 * crop%grass%zmow - mesh%dz(nodmow)
        do while (drz1 .gt. 0.d0)
          nodmow   = nodmow + 1
          drz1 = drz1 - mesh%dz(nodmow)
        enddo
      
      endif   
      
      if (swlossgrz .eq. 1) then
        
        ! Find node and layer for monitoring work-ability at start of grazing
        nodgrz = 1
        drz1       = -1.d0 * crop%grass%zgrz - mesh%dz(nodgrz)
        do while (drz1 .gt. 0.d0)
          nodgrz   = nodgrz + 1
          drz1 = drz1 - mesh%dz(nodgrz)
        enddo
         
      endif
      
      return

      case (2)

! === calculate potential rate and state variables ======================================

! --- rates of change of the grass variables ---------------------------------------------

      rid = dble(daycrop)
      
! --- check end of harvest
      if (crop%grass%flhrvendpot) then
        if (flearlyhrvendpot) then
          crop%grass%cropstartpot  = rid - 1.d0
        else
          crop%grass%cropstartpot  = rid
        endif
        crop%grass%pmowdm        = 0.d0
        crop%grass%pgrzdm        = 0.d0
        crop%wofost%plossdm       = 0.d0
      endif
      crop%grass%flhrvendpot      = .false.
      flearlyhrvendpot = .false.

! --- grass growth initiated by tsum from 1st day of calendar year
      crop%common%tsum = crop%common%tsum + max(0.0d0,atmo%Tav)
      if (.not. flGrassGrowth) then
        
        ! grass growth initiated by tsum
        if (crop%grass%swtsum.eq.1) then
          if (crop%common%tsum.ge.200.d0) then
            flGrassGrowth = .true.
          endif
        endif
        
        ! grass growth initiated by temperature, time and depth
        if (crop%grass%swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          if (dateGrassGrowth.eq.'undefined') call sumttd('dynamic',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif
      
        ! check if grass growth has started
        if (flGrassGrowth) then
          crop%grass%cropstartpot = rid
          crop%grass%cropstartact = rid
        endif

      endif
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.crop%grass%idregrpot) then

! ===   daily dry matter production ===

        gasspot = crop%wofost%pgasspot

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(crop%common%rmr*crop%wofost%wrtpot+crop%common%rml*crop%wofost%wlvpot+crop%common%rms*crop%wofost%wstpot)*afgen(crop%common%rfsetb,30,rid)
        teff = crop%common%q10**((atmo%Tav-25.0d0)/10.0d0)
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(crop%common%frtb,30,rid)
        fl = afgen(crop%common%fltb,30,rid)
        fs = afgen(crop%common%fstb,30,rid)
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
        cvf = 1.0d0/((fl/crop%common%cvl+fs/crop%common%cvs)*(1.0d0-fr)+fr/crop%common%cvr)
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
        if (crop%common%swrd.eq.3 .and. crop%wofost%wrtpot.gt.crop%common%wrtmax) then
          drrtpot = grrtpot
          drrtpot = max(drrtpot,crop%wofost%wrtpot*afgen (crop%common%rdrrtb,30,rid))
        else  
          drrtpot = crop%wofost%wrtpot*afgen (crop%common%rdrrtb,30,rid)
        endif  
        gwrtpot = grrtpot - drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/crop%kdif
        dslv2pot=crop%wofost%wlvpot*max(0.0d0,                                      &
     &                  min(0.03d0,0.03d0*(crop%common%laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = crop%common%ilvoldpot

        do while (restpot.gt.crop%common%lvpot(max(i1,1)).and.i1.ge.1)
          restpot = restpot-crop%common%lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (crop%common%lvagepot(max(i1,1)).gt.crop%common%span.and.restpot.gt.0.and.           &
     &                          i1.ge.1) then
          dalvpot = crop%common%lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.crop%common%lvagepot(max(i1,1)).gt.crop%common%span)
          dalvpot = dalvpot+crop%common%lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot   = dslvpot+dalvpot

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (crop%common%slatb,30,rid)
        if (crop%common%laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,atmo%Tav-crop%common%tbase)
          crop%common%glaiexpot = crop%common%laiexppot*crop%common%rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (crop%common%glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (crop%common%rdrstb,30,rid)*crop%wofost%wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        crop%grass%daygrowthpot = crop%grass%daygrowthpot + 1

!       Check trigger to start mowing event
        if (crop%grass%seqgrazmowpot(crop%grass%iseqgmpot) .eq. 2) then

          crop%grass%flHarvestpot = .false.
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (crop%wofost%tagppot .gt. dmharvest .or. (time%daynr .gt. daylastharvest  &
     &          .and. crop%wofost%tagppot .gt. dmlastharvest)) then
                crop%grass%flHarvestpot = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(crop%grass%dmmowtb,20,rid)
              if (crop%wofost%tagppot .gt. dmharvest .or.                           &
     &                 (crop%grass%daygrowthpot .gt. maxdaymow .and. crop%grass%iseqgmpot .gt. 1)) then
                crop%grass%flHarvestpot = .true.
              endif
            endif

          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(time%t1900 .gt. crop%grass%dateharvest(crop%grass%iharvest)) then
              crop%grass%flHarvestpot = .true.
            endif
          endif

!         In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
          if (crop%grass%flHarvestpot) then
            crop%grass%iseqgmpot = crop%grass%iseqgmpot + 1
            crop%common%slapot(1) = afgen (crop%common%slatb,30,rid)
            fl = afgen (crop%common%fltb,30,rid)
            fs = afgen (crop%common%fstb,30,rid)
            crop%wofost%wlvpot = crop%grass%mowrest / (1.d0 + (fs/fl))
            crop%wofost%wstpot = fs/fl*crop%wofost%wlvpot
            crop%wofost%dwlvpot = 0.0d0
            crop%wofost%dwstpot = 0.0d0
            crop%common%lvagepot(1) = 0.0d0
            crop%common%ilvoldpot = 1
            lasumpot = crop%wofost%wlvpot * crop%common%slapot(1)
            crop%common%laiexppot = lasumpot
            crop%common%lvpot(1) = crop%wofost%wlvpot
            
            gwstpot = 0.0d0
            gwrtpot = 0.0d0
            drlvpot = 0.0d0
            drstpot = 0.0d0
            drrtpot = 0.0d0
            
            crop%grass%daygrowthpot = 0
            
!           losses due to treading
            fralossmow = 0.d0
            if (swlossmow.eq.1) then
              fralossmow = afgen(crop%grass%lossmowtab,200,soil%h(nodmow))
            end if

!           harvest
            tagpspot = max(0.0d0,(crop%wofost%tagppot-(crop%wofost%wlvpot+crop%wofost%dwlvpot+crop%wofost%wstpot+crop%wofost%dwstpot)))
            crop%wofost%tagptpot = crop%wofost%tagptpot + tagpspot * (1.d0 - FraLossMow)

            crop%grass%cropendpot  = rid
            crop%grass%flhrvendpot = .true.
            crop%grass%pmowdm      = tagpspot * (1.d0 - FraLossMow)
            crop%wofost%plossdm     = tagpspot * FraLossMow
            
!           set regrowth delay
            idelaypot = int(afgen(crop%grass%DelayRegrowthTab,200,tagpspot))
            crop%grass%idregrpot = daycrop + idelaypot

          endif          
          
!       Check trigger to start grazing event          
        else if (crop%grass%seqgrazmowpot(crop%grass%iseqgmpot) .eq. 1 .or. crop%grass%seqgrazmowpot(crop%grass%iseqgmpot) .eq. 3) then

          crop%grass%flHarvestpot = .false.
          
          if (.not. crop%grass%flGrazingpot) then
              
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (crop%wofost%tagppot .gt. dmgrazing) then
                  crop%grass%flHarvestpot = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(crop%grass%dmgrztb,20,rid)
                if (crop%wofost%tagppot .gt. dmgrazing .or.                           &
     &                 (crop%grass%daygrowthpot .gt. maxdaygrz .and. crop%grass%iseqgmpot .gt. 1)) then
                  crop%grass%flHarvestpot = .true.
                endif
              endif
              
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(time%t1900 .gt. crop%grass%dateharvest(crop%grass%iharvest)) then
                crop%grass%flHarvestpot = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (crop%grass%flHarvestpot .or. crop%grass%flGrazingpot) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgrazpot = crop%grass%lsda(crop%grass%iseqgmpot) *                                 &
     &                         afgen(crop%grass%uptgrazingtab,200,crop%grass%lsda(crop%grass%iseqgmpot))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgrazpot = crop%grass%lsda(crop%grass%iseqgmpot) *                                &
     &                         afgen(crop%grass%lossgrazingtab,200,crop%grass%lsda(crop%grass%iseqgmpot))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(crop%grass%lossgrztab,200,soil%h(nodgrz))
            end if
            lossgrazpot = lossgrazpot + crop%wofost%tagppot * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. crop%grass%flGrazingpot) then
              crop%grass%daygrowthpot   = 0
              crop%grass%idaysgrazpot   = 0
              flDewoolingpot = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((crop%wofost%tagppot - uptgrazpot - lossgrazpot) .gt. tagprest) then
              
              crop%grass%flGrazingpot = .true.
              crop%common%cuptgrazpot  = crop%common%cuptgrazpot + uptgrazpot
          
!             distribute grazing over stems and leaves (living and dead parts)
              crop%wofost%wstpot  = crop%wofost%wstpot  - (uptgrazpot+lossgrazpot) * crop%wofost%wstpot  / crop%wofost%tagppot
              crop%wofost%dwstpot = crop%wofost%dwstpot - (uptgrazpot+lossgrazpot) * crop%wofost%dwstpot / crop%wofost%tagppot
              crop%wofost%dwlvpot = crop%wofost%dwlvpot - (uptgrazpot+lossgrazpot) * crop%wofost%dwlvpot / crop%wofost%tagppot
              grazlivinglvpot =   (uptgrazpot+lossgrazpot) * crop%wofost%wlvpot  / crop%wofost%tagppot
          
!             reduce leave weights
              i1 = crop%common%ilvoldpot
              do while (grazlivinglvpot .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglvpot .ge. crop%common%lvpot(i1)) then
                  grazlivinglvpot = grazlivinglvpot - crop%common%lvpot(i1)
                  crop%common%lvpot(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  crop%common%lvpot(i1) = crop%common%lvpot(i1) - grazlivinglvpot
                  grazlivinglvpot = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              crop%grass%cropendpot = rid
              crop%grass%pgrzdm     = crop%grass%pgrzdm + uptgrazpot
              crop%wofost%plossdm    = crop%wofost%tagppot * fralossgrz
              
!             Check number of days with grazing
              daysgrazpot  = int(afgen(crop%grass%daysgrazingtab,200,crop%grass%lsda(crop%grass%iseqgmpot)))
              crop%grass%idaysgrazpot = crop%grass%idaysgrazpot + 1
              if(crop%grass%idaysgrazpot .eq. daysgrazpot) then
                crop%grass%flGrazingpot = .false.
                crop%grass%flhrvendpot  = .true.
                if (crop%grass%seqgrazmowpot(crop%grass%iseqgmpot) .eq. 3) then
                  flDewoolingpot  = .true.
                endif
                crop%grass%daygrowthpot = 0
                crop%grass%iseqgmpot = crop%grass%iseqgmpot + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (crop%grass%flGrazingpot .or. swharvest .eq. 2) then
              crop%grass%flGrazingpot     = .false.
              crop%grass%flhrvendpot      = .true.
              flearlyhrvendpot = .true.
              if (crop%grass%seqgrazmowpot(crop%grass%iseqgmpot) .eq. 3 .and. crop%wofost%tagppot .gt. crop%grass%dewrest) then
                flDewoolingpot   = .true.
                flearlyhrvendpot = .false.
              endif
              crop%grass%daygrowthpot = 0
              crop%grass%iseqgmpot = crop%grass%iseqgmpot + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            crop%grass%idregrpot = daycrop

!           Dewooling after grazing event            
            if (flDewoolingpot) then

              crop%common%slapot(1) = afgen (crop%common%slatb,30,rid)
              fl = afgen (crop%common%fltb,30,rid)
              fs = afgen (crop%common%fstb,30,rid)
              crop%wofost%wlvpot = crop%grass%dewrest / (1.d0 + (fs/fl))
              crop%wofost%wstpot = fs/fl*crop%wofost%wlvpot
              crop%wofost%dwlvpot = 0.0d0
              crop%wofost%dwstpot = 0.0d0
              crop%common%lvagepot(1) = 0.0d0
              crop%common%ilvoldpot = 1
              lasumpot = crop%wofost%wlvpot * crop%common%slapot(1)
              crop%common%laiexppot = lasumpot
              crop%common%lvpot(1) = crop%wofost%wlvpot
              
              gwstpot = 0.0d0
              gwrtpot = 0.0d0
              drlvpot = 0.0d0
              drstpot = 0.0d0
              drrtpot = 0.0d0
              
!             Assumption: one day delay in regrowth after grazing
              crop%grass%idregrpot = daycrop + 1
              
            endif
            
          endif
          
        endif
        
        if (daycrop .ge. crop%grass%idregrpot) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(atmo%Tav-crop%common%tbase)/(35.0d0-crop%common%tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = crop%common%ilvoldpot
           do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.crop%common%lvpot(i1)) then
              dslvtpot = dslvtpot-crop%common%lvpot(i1)
              crop%common%lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              crop%common%lvpot(i1) = crop%common%lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (crop%common%lvagepot(max(i1,1)) .gt. crop%common%span .and. i1 .ge. 1)
              crop%common%lvpot(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          crop%common%ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = crop%common%ilvoldpot,1,-1
            crop%common%lvpot(i1+1) = crop%common%lvpot(i1)
            crop%common%slapot(i1+1) = crop%common%slapot(i1)
            crop%common%lvagepot(i1+1) = crop%common%lvagepot(i1)+fysdel*delt
          enddo
          crop%common%ilvoldpot = crop%common%ilvoldpot+1

! ---     new leaves in class 1
          crop%common%lvpot(1) = grlvpot*delt
          crop%common%slapot(1) = slatpot
          crop%common%lvagepot(1) = 0.d0

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          crop%wofost%wlvpot = 0.d0
          do i1 = 1,crop%common%ilvoldpot
            lasumpot = lasumpot+crop%common%lvpot(i1)*crop%common%slapot(i1)
            crop%wofost%wlvpot = crop%wofost%wlvpot+crop%common%lvpot(i1)
          enddo

          crop%common%laiexppot = crop%common%laiexppot+crop%common%glaiexpot*delt

        endif

! ---   dry weight of living plant organs
        crop%wofost%wrtpot = crop%wofost%wrtpot+gwrtpot*delt
        crop%wofost%wstpot = crop%wofost%wstpot+gwstpot*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        crop%wofost%dwrtpot = crop%wofost%dwrtpot+drrtpot*delt
        crop%wofost%dwlvpot = crop%wofost%dwlvpot+drlvpot*delt
        crop%wofost%dwstpot = crop%wofost%dwstpot+drstpot*delt

! ---   dry weight of dead and living plant organs
        twlvpot = crop%wofost%wlvpot+crop%wofost%dwlvpot
        twstpot = crop%wofost%wstpot+crop%wofost%dwstpot
        crop%wofost%tagppot = twlvpot+twstpot

! ---   leaf area index
        crop%common%laipot = lasumpot+crop%common%ssa*crop%wofost%wstpot
!       prevent immediate lai reduction at emergence
!       KRO-BOO-20160403: suppressed because deviates from Wofost
!       laipot = max(laipot, laiem)

        ! root extension
        if (crop%common%swrd.eq.1) then
          crop%common%rdpot = afgen (crop%common%rdtb,22,rid)
          crop%common%rdpot = min(crop%common%rdpot,crop%common%rdm)
        elseif (crop%common%swrd.eq.2) then
          rrpot = min (crop%common%rdm-crop%common%rdpot,crop%common%rri)
          if (fr.le.0.0d0 .or. crop%wofost%pgasspot.lt.1.0d0) rrpot = 0.0d0
          crop%common%rdpot = crop%common%rdpot + rrpot
        elseif (crop%common%swrd.eq.3) then
          crop%common%rdpot = afgen (crop%common%rlwtb,22,crop%wofost%wrtpot)
          crop%common%rdpot = min(crop%common%rdpot,crop%common%rdm)
        endif

      endif

      ! mirror grass case(2) potential state
      crop%common%cuptgrazpot   = crop%common%cuptgrazpot

      return

      case (3)

! === calculate actual rate and state variables ======================================

! --- check end of harvest
      if (crop%grass%flhrvendact) then
        if (flearlyhrvendact) then
          crop%grass%cropstartact  = rid - 1.d0
        else
          crop%grass%cropstartact  = rid
        endif
        crop%grass%mowdm        = 0.d0
        crop%grass%grzdm        = 0.d0
        crop%wofost%lossdm       = 0.d0
      endif
      crop%grass%flhrvendact      = .false.
      flearlyhrvendact = .false.

! --- rates of change of the crop variables ---------------------------------------------
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop .ge. crop%grass%idregr) then

! ===   daily dry matter production ===

! ---   water stress reduction of pgass to gass
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if(dabs(atmo%ptra).lt.nihil) then
          crop%common%reltr = 1.0d0
        else
          crop%common%reltr = max(0.0d0,min(1.0d0,soil%tra/atmo%ptra))
        endif
        gass = crop%wofost%pgass * crop%common%reltr

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (crop%common%rmr*crop%wofost%wrt+crop%common%rml*crop%wofost%wlv+crop%common%rms*crop%wofost%wst)*afgen(crop%common%rfsetb,30,rid)
        teff = crop%common%q10**((atmo%Tav-25.0d0)/10.0d0)
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   partitioning factors (relevant for restart)
        fr = afgen(crop%common%frtb,30,rid)
        fl = afgen(crop%common%fltb,30,rid)
        fs = afgen(crop%common%fstb,30,rid)

! ---   dry matter increase
        cvf = 1.0d0/((fl/crop%common%cvl+fs/crop%common%cvs)*(1.0d0-fr)+fr/crop%common%cvr)
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
        if (crop%common%swrd.eq.3 .and. soil%flWrtNonox) grrt = 0.d0 
        if (crop%common%swrd.eq.3 .and. crop%wofost%wrt.gt.crop%common%wrtmax) then
          drrt = grrt
          drrt = max(drrt,crop%wofost%wrt*afgen (crop%common%rdrrtb,30,rid))
        else  
          drrt = crop%wofost%wrt*afgen (crop%common%rdrrtb,30,rid)
        endif  
        gwrt = grrt-drrt

! ---   growth rate leaves

! ---   weight of new leaves
        admi = (1.0d0-fr)*dmi        
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        dslv1 = crop%wofost%wlv*(1.0d0-crop%common%reltr)*crop%common%perdl
        laicr = 3.2d0/crop%kdif
        dslv2 = crop%wofost%wlv*max(0.0d0,min(0.03d0,0.03d0*(crop%lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = crop%common%ilvold

        do while (rest.gt.crop%common%lv(max(i1,1)).and.i1.ge.1)
          rest = rest-crop%common%lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (crop%common%lvage(max(i1,1)).gt.crop%common%span.and.rest.gt.0.and.i1.ge.1) then
          dalv = crop%common%lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.crop%common%lvage(max(i1,1)).gt.crop%common%span)
          dalv = dalv+crop%common%lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv   = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (crop%common%slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (crop%common%laiexp.lt.6.0d0) then
          dteff = max (0.0d0,atmo%Tav-crop%common%tbase)
          crop%common%glaiex = crop%common%laiexp*crop%common%rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (crop%common%glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        drst1 = crop%wofost%wst*(1.0d0-crop%common%reltr)*crop%common%perdl
! ---   death of stems due to ageing
        drst2 = afgen (crop%common%rdrstb,30,rid)*crop%wofost%wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        crop%grass%daygrowth = crop%grass%daygrowth + 1

!       Check trigger to start mowing event
        if (crop%grass%seqgrazmow(crop%grass%iseqgm) .eq. 2) then
          
          crop%grass%flHarvest = .false.   
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (crop%wofost%tagp .gt. dmharvest .or. (time%daynr .gt. daylastharvest  &
     &          .and. crop%wofost%tagp .gt. dmlastharvest)) then
                crop%grass%flHarvest = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(crop%grass%dmmowtb,20,rid)
              if (crop%wofost%tagp .gt. dmharvest .or.                           &
     &                 (crop%grass%daygrowth .gt. maxdaymow .and. crop%grass%iseqgm .gt. 1)) then
                crop%grass%flHarvest = .true.
              endif
            endif
          
          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(time%t1900 .gt. crop%grass%dateharvest(crop%grass%iharvest)) then
              crop%grass%iharvest = crop%grass%iharvest + 1
              crop%grass%flHarvest = .true.
            endif
          endif
          
!       In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
        if (crop%grass%flHarvest) then
          crop%grass%iseqgm = crop%grass%iseqgm + 1
          crop%common%sla(1) = afgen (crop%common%slatb,30,rid)
          fl = afgen (crop%common%fltb,30,rid)
          fs = afgen (crop%common%fstb,30,rid)
          crop%wofost%wlv = crop%grass%mowrest / (1.d0 + (fs/fl))
          crop%wofost%wst = fs/fl*crop%wofost%wlv
          crop%wofost%dwlv = 0.0d0
          crop%wofost%dwst = 0.0d0
          crop%common%lvage(1) = 0.0d0
          crop%common%ilvold = 1
          lasum = crop%wofost%wlv * crop%common%sla(1)
          crop%common%laiexp = lasum
          crop%common%lv(1) = crop%wofost%wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          crop%grass%daygrowth = 0

!         losses due to treading
          FraLossMow = 0.d0
          if (swlossmow.eq.1) then
            FraLossMow = afgen(crop%grass%lossmowtab,200,soil%h(nodmow))
          end if
          
!         harvest
          tagps = max (0.0d0,(crop%wofost%tagp-(crop%wofost%wlv+crop%wofost%dwlv+crop%wofost%wst+crop%wofost%dwst)))
          crop%wofost%tagpt = crop%wofost%tagpt + tagps * (1.d0 - fralossmow)

          crop%grass%cropendact  = rid
          crop%grass%flhrvendact = .true.
          crop%grass%mowdm   = tagps * (1.d0 - FraLossMow)
          crop%wofost%lossdm  = tagps * FraLossMow
          
! ---     set regrowth delay
          idelay = int(afgen(crop%grass%DelayRegrowthTab,200,tagps))
          crop%grass%idregr = daycrop + idelay

        endif
          
!       Check trigger to start grazing event          
        else if (crop%grass%seqgrazmow(crop%grass%iseqgm) .eq. 1 .or. crop%grass%seqgrazmow(crop%grass%iseqgm) .eq. 3) then

          crop%grass%flHarvest = .false.
            
          if (.not. crop%grass%flGrazing) then
            
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (crop%wofost%tagp .gt. dmgrazing) then
                  crop%grass%flHarvest = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(crop%grass%dmgrztb,20,rid)
                if (crop%wofost%tagp .gt. dmgrazing .or.                           &
     &            (crop%grass%daygrowth .gt. maxdaygrz .and. crop%grass%iseqgm .gt. 1)) then
                  crop%grass%flHarvest = .true.
                endif
              endif
            
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(time%t1900 .gt. crop%grass%dateharvest(crop%grass%iharvest)) then
                crop%grass%iharvest = crop%grass%iharvest + 1
                crop%grass%flHarvest = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (crop%grass%flHarvest .or. crop%grass%flGrazing) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = crop%grass%lsda(crop%grass%iseqgm)*afgen(crop%grass%uptgrazingtab,200,crop%grass%lsda(crop%grass%iseqgm))            

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = crop%grass%lsda(crop%grass%iseqgm) *                                &
     &                         afgen(crop%grass%lossgrazingtab,200,crop%grass%lsda(crop%grass%iseqgm))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(crop%grass%lossgrztab,200,soil%h(nodgrz))
            end if
            lossgraz = lossgraz + crop%wofost%tagp * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. crop%grass%flGrazing) then
              crop%grass%daygrowth   = 0
              crop%grass%idaysgraz   = 0
              flDewooling = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((crop%wofost%tagp - uptgraz - lossgraz) .gt. tagprest) then
              
              crop%grass%flGrazing = .true.
              crop%common%cuptgraz  = crop%common%cuptgraz + uptgraz

!             distribute grazing over stems and leaves (living and dead parts)
              crop%wofost%wst  = crop%wofost%wst  -  (uptgraz+lossgraz) * crop%wofost%wst  / crop%wofost%tagp
              crop%wofost%dwst = crop%wofost%dwst -  (uptgraz+lossgraz) * crop%wofost%dwst / crop%wofost%tagp
              crop%wofost%dwlv = crop%wofost%dwlv -  (uptgraz+lossgraz) * crop%wofost%dwlv / crop%wofost%tagp
              grazlivinglv = (uptgraz+lossgraz) * crop%wofost%wlv  / crop%wofost%tagp
          
!             reduce leave weights
              i1 = crop%common%ilvold
              do while (grazlivinglv .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglv .ge. crop%common%lv(i1)) then
                  grazlivinglv = grazlivinglv - crop%common%lv(i1)
                  crop%common%lv(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  crop%common%lv(i1) = crop%common%lv(i1) - grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              crop%grass%cropendact = rid
              crop%grass%grzdm      = crop%grass%grzdm + uptgraz
              crop%wofost%lossdm     = crop%wofost%tagp * fralossgrz
              
!             Check number of days with grazing
              daysgraz  = int(afgen(crop%grass%daysgrazingtab,200,crop%grass%lsda(crop%grass%iseqgm)))
              crop%grass%idaysgraz = crop%grass%idaysgraz + 1
              if(crop%grass%idaysgraz .eq. daysgraz) then
                crop%grass%flGrazing   = .false.
                crop%grass%flhrvendact = .true.
                if (crop%grass%seqgrazmow(crop%grass%iseqgm) .eq. 3) then
                  flDewooling  = .true.
                endif
                crop%grass%daygrowth = 0
                crop%grass%iseqgm = crop%grass%iseqgm + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (crop%grass%flGrazing .or. swharvest .eq. 2) then
              crop%grass%flGrazing        = .false.
              crop%grass%flhrvendact      = .true.
              flearlyhrvendact = .true.
              if (crop%grass%seqgrazmow(crop%grass%iseqgm) .eq. 3 .and. crop%wofost%tagp .gt. crop%grass%dewrest) then
                flDewooling      = .true.
                flearlyhrvendact = .false.
              endif
              crop%grass%daygrowth = 0
              crop%grass%iseqgm = crop%grass%iseqgm + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            crop%grass%idregr = daycrop

!           Dewooling after grazing event            
            if (flDewooling) then

              crop%common%sla(1) = afgen (crop%common%slatb,30,rid)
              fl = afgen (crop%common%fltb,30,rid)
              fs = afgen (crop%common%fstb,30,rid)
              crop%wofost%wlv = crop%grass%dewrest / (1.d0 + (fs/fl))
              crop%wofost%wst = fs/fl*crop%wofost%wlv
              crop%wofost%dwlv = 0.0d0
              crop%wofost%dwst = 0.0d0
              crop%common%lvage(1) = 0.0d0
              crop%common%ilvold = 1
              lasum = crop%wofost%wlv * crop%common%sla(1)
              crop%common%laiexp = lasum
              crop%common%lv(1) = crop%wofost%wlv
    
              gwst = 0.0d0
              gwrt = 0.0d0
              drlv = 0.0d0
              drst = 0.0d0
              drrt = 0.0d0
    
!             Assumption: one day delay in regrowth after grazing
              crop%grass%idregr = daycrop + 1

            endif
            
          endif
          
        endif

        if (daycrop .ge. crop%grass%idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(atmo%Tav-crop%common%tbase)/(35.0d0-crop%common%tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = crop%common%ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.crop%common%lv(i1)) then
              dslvt = dslvt-crop%common%lv(i1)
              crop%common%lv(i1) = 0.0d0
              i1 = i1-1
            else
              crop%common%lv(i1) = crop%common%lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (crop%common%lvage(max(i1,1)).gt.crop%common%span.and.i1.ge.1)
              crop%common%lv(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          crop%common%ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = crop%common%ilvold,1,-1
            crop%common%lv(i1+1) = crop%common%lv(i1)
            crop%common%sla(i1+1) = crop%common%sla(i1)
            crop%common%lvage(i1+1) = crop%common%lvage(i1)+fysdel*delt
          enddo
          crop%common%ilvold = crop%common%ilvold+1

! ---     new leaves in class 1
          crop%common%lv(1) = grlv*delt
          crop%common%sla(1) = slat
          crop%common%lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          crop%wofost%wlv = 0.d0
          do i1 = 1,crop%common%ilvold
            lasum = lasum+crop%common%lv(i1)*crop%common%sla(i1)
            crop%wofost%wlv = crop%wofost%wlv+crop%common%lv(i1)
          enddo

          crop%common%laiexp = crop%common%laiexp+crop%common%glaiex*delt

        endif

! ---   dry weight of living plant organs
        crop%wofost%wrt = crop%wofost%wrt+gwrt*delt
        crop%wofost%wst = crop%wofost%wst+gwst*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        crop%wofost%dwrt = crop%wofost%dwrt+drrt*delt
        crop%wofost%dwlv = crop%wofost%dwlv+drlv*delt
        crop%wofost%dwst = crop%wofost%dwst+drst*delt

! ---   dry weight of dead and living plant organs
        twlv = crop%wofost%wlv+crop%wofost%dwlv
        twst = crop%wofost%wst+crop%wofost%dwst
        crop%wofost%tagp = twlv+twst

! ---   leaf area index
        crop%lai = lasum+crop%common%ssa*crop%wofost%wst
        crop%common%laimax = max (crop%lai,crop%common%laimax)

! ---   update normalized cumulative root density based on root extraction or stress (cumdens)
        if (crop%common%swrdc .eq. 1) call update_rootdistribution(state)
        
        ! root extension
        if (crop%common%swrd.eq.1) then
          crop%common%rd = afgen (crop%common%rdtb,22,rid)
          crop%common%rd = min(crop%common%rd,crop%common%rdm)
        elseif (crop%common%swrd.eq.2) then
          rr = min (crop%common%rdm-crop%common%rd,crop%common%rri)
          if (fr.le.0.0d0 .or. crop%wofost%pgass.lt.1.0d0 .or.                    &
     &        soil%flWrtNonox) rr = 0.0d0 
          if (crop%common%swdmi2rd.eq.1 .and. crop%wofost%pgass.ge.1.0d0)              rr = rr * gass/crop%wofost%pgass
          crop%common%rd = crop%common%rd + rr
        elseif (crop%common%swrd.eq.3) then
          crop%common%rd = afgen (crop%common%rlwtb,22,crop%wofost%wrt)
          crop%common%rd = min(crop%common%rd,crop%common%rdm)
        endif

! ---   set crop height and cropfactor
        if (crop%swcf.ne.3) then
          crop%common%cf = afgen (crop%fixed%cftb,(2*magrs),rid)
          crop%common%ch = afgen (crop%fixed%chtb,(2*magrs),rid)
        else
          crop%common%cf = afgen (crop%fixed%cftb,(2*magrs),crop%lai)
          crop%fixed%cfeic = afgen (crop%fixed%cfeictb,(2*magrs),crop%lai)
          crop%common%ch = afgen(crop%fixed%chtb,(2*magrs),crop%lai)
        endif
        ! cfeic write retired (was: crop%fixed%cfeic = cfeic)

! ---   update canopy storage capacity
        if (crop%common%swinter.eq.3) then
          atmo%siccapact = 0.0d0   ! [GR-CROP 2026-05-25] swinter=3 stub-errored on TOML; siccaplai always 0
        endif

      endif

      ! mirror grass case(3) actual state
      crop%wofost%tagp        = crop%wofost%tagp
      crop%common%cuptgraz    = crop%common%cuptgraz

      return

      case default
         call fatalerr_collected ('Grass', 'Illegal value for TASK')
      end select

      end associate
      return
      end subroutine grass

      end module cropgrass_runtime_mod
