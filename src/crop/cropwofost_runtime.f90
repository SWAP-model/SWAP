! cropwofost_runtime.f90 — type-2 (WOFOST) crop runtime dispatcher.
!
! [GR-CROP 2026-05-25] crop-sweep:
!   - macp/magrs sourced from swap_array_dimensions.
!   - rdmax read via crop_cfg%rdmax (Class B); daycrop read via crop%common%daycrop (Class A).
!   - swbulb reads → state%crop%wofost%swbulb (logical mirror; default false; cfg.bulb=1 stub-errored on TOML).
!   - Wofost working state (gasst/gasstpot/mrest/mrestpot/tadw/tadwpot/fbl/drbl/drblpot)
!     migrated from legacy globals to local SAVE (module-level `save` carries state across task=1..4).
!   - 12 single-file config snapshots (rdrns/dvsnlt/dvsnt/fntrt/tcnt + fraharlosorm_lv/st/so +
!     vernbase/verndvs/vernsat/vernrtb) migrated to public module SAVE (cw_*) in cropwofost_init_mod.
!   - Remaining legacy globals: cross-file with cropgrowth dispatcher (Task 9) or other
!     crop sub-arcs (cropgrass/cropfixed/oxygenstress/rootextraction); retired there.
! ----------------------------------------------------------------------
      module cropwofost_runtime_mod
      implicit none
      private

      public :: wofost
      public :: wofost_apply_nstress

      contains

! ----------------------------------------------------------------------
      subroutine wofost_apply_nstress(state)
! ----------------------------------------------------------------------
! Purpose: Apply WOFOST nitrogen-stress reduction to potential gross
!          assimilation. Called by cropgrowth.f90 dispatcher when
!          state%crop%common%flCropNut is true. The nutrient cluster
!          (anlv/anst/nmxlv/nmaxlv/.../nlue/lrnr/lsnr/fstr/...) still
!          lives in legacy_state — this wrapper isolates the bare-global
!          access so the dispatcher can stay `use variables`-free.
! ----------------------------------------------------------------------
      ! [GR-CROP 2026-05-25] Nutrient cluster now lives as module-level cw_* SAVE
      ! in cropwofost_init_mod (single-file scope across the cropwofost pair).
      use cropwofost_init_mod, only: cw_nlue, cw_anlv, cw_anst, cw_nmxlv, &
                                     cw_nmaxlv, cw_nmaxst, cw_nmaxrt,     &
                                     cw_lrnr, cw_lsnr, cw_nni, cw_rnflv,  &
                                     cw_rnfst, cw_frnx, cw_fstr
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      call NUTRIE(cw_NLUE, state%crop%wofost%wlv, state%crop%wofost%wst, &
                  state%crop%common%dvs, cw_ANLV, cw_ANST, cw_NMXLV, cw_NMAXLV,   &
                  cw_NMAXST, cw_NMAXRT, cw_LRNR, cw_LSNR, cw_NNI, cw_RNFLV, cw_RNFST, cw_FRNX, cw_FSTR)
      state%crop%wofost%pgass = state%crop%wofost%pgass * cw_FSTR

      return
      end subroutine wofost_apply_nstress

      subroutine wofost(task, state)
! ----------------------------------------------------------------------
!     update             : march 2015
!     date               : october 2004
!     purpose            : detailed crop growth routine
! ----------------------------------------------------------------------
      use swap_log, only: log_warn
      use swap_array_dimensions, only: macp, magrs
      ! [GR-CROP 2026-05-25] crop-sweep:
      !   - macp/magrs sourced from swap_array_dimensions.
      !   - rdmax read via crop_cfg%rdmax (Class B direct read).
      !   - daycrop read via crop%common%daycrop (Class A state-rebind).
      !   - Remaining legacy globals: writer/reader targets feeding cropgrowth
      !     dispatcher (Task 9) and cross-file consumers (cropgrass/cropfixed/
      !     oxygenstress/rootextraction); cannot retire here.
      ! [GR-CROP 2026-05-25] outfil/pathwork/project read via state%cfg%general
      ! at each outbalcrop* call site; file is now `use variables`-free.
      use wofost_soil_interface
      ! [GR-CROP 2026-05-25] cw_* snapshots — written by cropwofost_init_mod%apply_cropwofost_nutrient
      ! Nutrient cluster + harvest/vernalisation fractions. Single-file scope.
      use cropwofost_init_mod, only: cw_rdrns, cw_dvsnlt, cw_dvsnt, cw_fntrt, cw_tcnt,            &
                                     cw_fraharlosorm_lv, cw_fraharlosorm_so, cw_fraharlosorm_st, &
                                     cw_vernbase, cw_verndvs, cw_vernsat, cw_vernrtb,            &
                                     cw_nlue, cw_lrnr, cw_lsnr, cw_rnflv, cw_rnfst, cw_rnfrt,     &
                                     cw_frnx, cw_nlai, cw_nmaxso, cw_npart, cw_nfixf, cw_nsla,    &
                                     cw_nmxlv, cw_ilnmxl,                                         &
                                     cw_anlv, cw_anst, cw_nni, cw_fstr,                          &
                                     cw_nmaxlv, cw_nmaxst, cw_nmaxrt
      use array_utils, only: interpol, afgen, insw
      use soilhydraulics_utils, only: watcon
        use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      ! GR-CROPWS Phase 0: update_rootdistribution extracted to cropgrowth_helpers_mod
      use cropgrowth_helpers_mod, only: update_rootdistribution
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer   i1,task,swhydrlift,i

      real(8)   asrc,ccheck,cvf
      real(8)   laicr,lasum,mres
      ! [GR-CROP 2026-05-25] wofost working state migrated from legacy globals to local SAVE.
      ! These are single-file persistent state across task=1/2/3/4 calls — module 'save' below.
      real(8)   gasst, gasstpot, mrest, mrestpot
      real(8)   tadw, tadwpot
      real(8)   fbl, drbl, drblpot
      real(8)   dalv,delt,dmi
      real(8)   drrt,drst,dslv,dteff,dtsum,dvr
      real(8)   dvred,fl,fo,fr,fs,drlv
      real(8)   fysdel,gass,gla,grlv,grrt,grst,grso
      real(8)   gwso,gwst,rmres
      real(8)   slat,teff,twlv,twst
      real(8)   lasumpot
      real(8)   gasspot,rmrespot,mrespot,asrcpot,dmipot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot
      real(8)   dslvpot,restpot,dalvpot,drlvpot,gwsopot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot
      real(8)   dummy
      character(len=200) messag
      
! --- rooting
      real(8)   rrpot,rr
      
! --- only for bulb crops (tulips etc..)
      real(8)   admi,grbl,gwbl,grblpot,gwblpot
      real(8)   twbl,twblpot

      real(8) Fstress
      integer nut

      real(8) NDEMTO
      real(8) ANRT,ANSO,ANLVI,ANSTI,ANRTI,ANSOI !,cw_ANLV,cw_ANST
      real(8) RNLDRT,RNLDST,RNLDLV,NDEML,NDEMS,NDEMR,NDEMSO
      real(8) NSUPSO,ATN,RNSO,NLIMIT,NUPTR
      real(8) NFIXTR,ATNLV,ATNST,ATNRT,RNTLV,RNTST,RNTRT
      real(8) RNULV,RNUST,RNURT,NUPTT,NFIXTT
      real(8) RNLV,RNST,RNRT
      real(8) NLOSSL,NLOSSR,NLOSSS  
      real(8) NBALAN

      real(8) NdemandBioFix
      real(8) ombalan,wlvt0,wstt0,wsot0,wrtt0,storagediff,drso
      real(8) FraDeceasedLvToSoil
      real(8) HarLosOrm_rt, HarLosOrm_lv, HarLosOrm_st, HarLosOrm_so
      real(8) HarLosNit_rt, HarLosNit_lv, HarLosNit_st, HarLosNit_so
      real(8) HarLosOrm_dwst, HarLosOrm_dwso, HarLosOrm_dwlv

      real(8) HarLosNit_dwst, HarLosNit_dwso, HarLosNit_dwlv
      real(8) idwlvCrop,idwlvSoil, NLOSSLDeceasedLvToSoil
      character(len=200) filnam

! --- only for soybean
      integer swsoybean
      logical flrfphotoveg
      real(8) mg,dvsi,dvrmax1,dvrmax2,tmaxdvr,tmindvr,toptdvr
      real(8) rfmgphotop, rfmgtemp
      logical flphenodayl
      real(8) popt
      real(8) pcrt

! --- only for bulb crops (tulips etc..)
      real(8) wstem

      parameter (delt=1.0d0)

! --- only for vernalisation
      logical flvernalised
      real(8) r,vern,vernfac,vernrate

      save
! ----------------------------------------------------------------------
      ! [GR-CROP 2026-05-25] sub-record associate (crop-sweep convention).
      associate( crop     => state%crop,            &
                 soil     => state%soilwater,       &
                 mesh     => state%mesh,            &
                 atmo     => state%atmosphere,      &
                 time     => state%timecontrol,     &
                 crop_cfg => state%cfg%crop         )

      select case (task)

      case (1)

! === initialization ====================================================

! --- read general crop data: dispatch on per-rotation typed-config cache
!     (ADR 0016). Falls back to legacy reader for rotations whose
!     .crp.toml is not yet authored. Teardown: end of Phase 4 removes
!     the else-branch.
      block
         use cropwofost_init_mod, only: cropwofost_init_from_config
         logical :: use_cache
         use_cache = .false.
         if (allocated(crop_cfg%rotation_loaded)) then
            if (crop%common%icrop >= 1 .and. crop%common%icrop <= size(crop_cfg%rotation_loaded)) then
               if (crop_cfg%rotation_loaded(crop%common%icrop)) use_cache = .true.
            end if
         end if
         if (use_cache) then
            call cropwofost_init_from_config(crop_cfg%rotation_wofost(crop%common%icrop), &
                                             crop%common%icrop, FraDeceasedLvToSoil, state)
            ! swhydrlift is read by legacy readwofost only inside swdrought=2
            ! branch (stub-errored in Phase 2). Set to 0 here to mirror the
            ! default.
            swhydrlift = 0
         else
            ! ADR 0016 cache-miss: typed config required for type=2 rotations.
            ! No silent legacy fallback — the user must author cropwofost.crp.toml.
            call fatalerr_collected('cropgrowth/Wofost', &
               'cropwofost rotation has no loaded .crp.toml — author the file or use the legacy executable.')
         end if
      end block

! --- if crop based on calendar is still active, but already harvested
      if (crop%common%flCropHarvest) return

! --- n-p-k
      if( crop%common%flCropNut) then
!        Legacy nutrient parameters (cw_LRNR, cw_LSNR, cw_NLAI, cw_NLUE, cw_NMAXSO,
!        cw_NPART, cw_NFIXF, cw_NSLA, cw_RNFLV/RT/ST, cw_tcnt, cw_dvsnlt, cw_dvsnt, RDRNS,
!        cw_fntrt, cw_FRNX, cw_NMXLV, cw_fraharlosorm_lv/st/so) used to be read here
!        from <cropfil>.crp via TTutil rdinit/rdsdou. Read block deleted
!        as part of legacy readers physical deletion. These globals are
!        now populated by apply_cropwofost_nutrient when crop%common%flCropNut=true
!        on the active rotation (ADR 0025 N1, ADR 0028 N3).

!        open output files and write header
         if (crop%common%icrop.eq.1) then
            call outbalcropOM1(1,state%cfg%general%pathwork,state%cfg%general%outfil,state%cfg%general%project,time%date,crop%common%daycrop,  &
     &         time%t,crop%common%dvs,crop%common%tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)
            call outbalcropOM2(1,state%cfg%general%pathwork,state%cfg%general%outfil,state%cfg%general%project,time%date,crop%common%daycrop,  &
     &         time%t,crop%common%dvs,crop%common%tsum,storagediff,crop%wofost%wlv,crop%wofost%wst,crop%wofost%wso,crop%wofost%wrt,delt,         &
     &         grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)
            call outbalcropN(1,state%cfg%general%pathwork,state%cfg%general%outfil,state%cfg%general%project,time%date,crop%common%daycrop, &
     &         time%t,crop%common%dvs,crop%common%tsum,nuptt,nfixtt,anlvi,ansti,anrti,ansoi,cw_anlv,&
     &         cw_anst,anrt,anso,nlossl,nlossr,nlosss,nbalan,cw_nni)
         endif
      endif

! --- maximum rooting depth
      if (crop%common%swrd.eq.1) then
        crop%common%rdm = crop_cfg%rdmax
      elseif (crop%common%swrd.eq.2) then
        crop%common%rdm = min(crop_cfg%rdmax,crop%common%rdc)
      elseif (crop%common%swrd.eq.3) then
        crop%common%rdc = afgen (crop%common%rlwtb,22,crop%common%wrtmax)
        crop%common%rdm = min(crop_cfg%rdmax,crop%common%rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (time%t1900 - time%tstart .gt. tiny .or. soil%swinco .ne. 3 .or.           &
     &   dabs(time%t1900 - crop%common%cropstart) .lt. tiny) then

        crop%common%dvs = 0.0d0
        crop%wofost%flanthesis = .false.
        crop%common%tsum = 0.0d0
        fr = afgen (crop%common%frtb,30,crop%common%dvs)
        fl = afgen (crop%common%fltb,30,crop%common%dvs)
        fs = afgen (crop%common%fstb,30,crop%common%dvs)
        fo = afgen (crop%common%fotb,30,crop%common%dvs)
! --- only for bulb crops (tulips etc..)
        if(crop%wofost%swbulb) then
           fbl = afgen (crop%common%fbltb,30,crop%common%dvs)
           crop%wofost%plwt = crop%wofost%plwti
        endif
        crop%common%sla(1) = afgen (crop%common%slatb,30,crop%common%dvs)
        crop%common%lvage(1) = 0.0d0
        crop%common%ilvold = 1
        crop%common%slapot(1) = afgen (crop%common%slatb,30,crop%common%dvs)
        crop%common%lvagepot(1) = 0.0d0
        crop%common%ilvoldpot = 1

! ---   initial state variables of the crop
        crop%wofost%wrt = fr*crop%common%tdwi
        crop%wofost%wrtmin = crop%wofost%wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        crop%wofost%wrtpot = crop%wofost%wrt
        tadw = (1.0d0-fr)*crop%common%tdwi
        tadwpot = tadw
        crop%wofost%wst = fs*tadw
        crop%wofost%wstpot = crop%wofost%wst
        crop%wofost%wso = fo*tadw
        crop%wofost%wsopot = crop%wofost%wso
        crop%wofost%wlv = fl*tadw
        crop%wofost%wlvpot = crop%wofost%wlv
! --- only for bulb crops (tulips etc..)
        if(crop%wofost%swbulb) then
!          blad bij opkomst is ondergronds: lai vanuit ingelezen laiem,
!          crop%common%sla(l) aangepast aan initieel bladgewicht en laiem,
!          stengelgewicht bij opkomst niet meegenomen bij lai-berekening
           crop%common%sla(1) = crop%common%laiem / crop%wofost%wlv
           wstem = crop%wofost%wst
           crop%wofost%wbl = fbl*tadw
           crop%wofost%wblpot = crop%wofost%wbl
        else
!          KRO-BOO-20160403: intro because comparison with Wofost
           crop%common%laiem = crop%wofost%wlv*crop%common%sla(1)
        endif
        crop%common%lv(1) = crop%wofost%wlv
        crop%common%lvpot(1) = crop%wofost%wlv
        lasum = crop%common%laiem
        lasumpot = crop%common%laiem
        crop%common%laiexp = crop%common%laiem
        crop%common%laiexppot = crop%common%laiem
        crop%common%glaiex = 0.0d0
        crop%common%glaiexpot = 0.0d0
        crop%common%laimax = crop%common%laiem
! --- only for bulb crops (tulips etc..)
        if(crop%wofost%swbulb) then
            crop%lai = lasum+crop%common%ssa*(crop%wofost%wst-wstem)+crop%common%spa*crop%wofost%wso
            crop%wofost%dwbl = 0.0d0
            crop%wofost%dwblpot = 0.0d0
        else
            crop%lai = lasum+crop%common%ssa*crop%wofost%wst+crop%common%spa*crop%wofost%wso
        endif
        crop%common%laipot = crop%lai
        crop%wofost%dwrt = 0.0d0
        crop%wofost%dwrtpot = 0.0d0
        crop%wofost%dwlv = 0.0d0
        crop%wofost%dwlvCrop = 0.0d0
        crop%wofost%dwlvSoil = 0.0d0
        crop%wofost%dwlvpot = 0.0d0
        crop%wofost%dwso = 0.0d0
        crop%wofost%dwst = 0.0d0
        crop%wofost%dwstpot = 0.0d0
        if(crop%common%flCropNut) then
          WLVt0 = crop%wofost%wlv
          WSTt0 = crop%wofost%wst
          WSOt0 = crop%wofost%wso
          WRTt0 = crop%wofost%wrt
        endif
        
! ---   n-p-k 
        if( crop%common%flCropNut) then
!******************************************************************
!         initial maximum nutrient concentrations in plant organs 
!         per kg biomass [kg N kg-1 dry biomass] at sowing added IS
!******************************************************************        
          call nutrsow(cw_anlv,cw_anst,anrt,anso)
!******************************************************************
!         initial maximum nutrient concentrations in plant organs 
!         per kg biomass [kg N kg-1 dry biomass] at emergence added IS
!******************************************************************
          call nutremrg(cw_nmxlv,cw_lsnr,cw_lrnr,crop%wofost%wlv,crop%wofost%wst,crop%wofost%wrt,                    &
     &      cw_anlv,cw_anst,anrt,anso,anlvi,ansti,anrti,ansoi,crop%common%dvs)
        endif

! --- actual rooting depth
        if (crop%common%swrd.eq.1) then
          crop%common%rd = afgen (crop%common%rdtb,22,crop%common%dvs)
          crop%common%rd = min(crop%common%rd,crop%common%rdm)
        elseif (crop%common%swrd.eq.2) then
          crop%common%rd = min(crop%common%rdi,crop%common%rdm)
        elseif (crop%common%swrd.eq.3) then
          crop%common%rdi = afgen (crop%common%rlwtb,22,crop%wofost%wrt)
          crop%common%rd = min(crop%common%rdi,crop%common%rdm)
        endif
        crop%common%rdpot = crop%common%rd
        
! --- initial summation variables of the crop
        gasst = 0.0d0
        gasstpot = 0.0d0
        mrest = 0.0d0 
        mrestpot = 0.0d0 
        crop%wofost%cwdm = 0.0d0
        crop%wofost%cwdmpot = 0.0d0
! --- only for vernalisation
        vern = 0.0d0             ! vernalisation state (d)
        flvernalised = .FALSE.   ! crop not vernalised (-)

        ! [GR-CROP 2026-05-25] swbulb legacy-global mirror dropped — state%crop%wofost%swbulb
        ! is already seeded false by swap_state init; cfg%bulb%swbulb=1 is stub-errored on TOML.

! --- end skip above initialization if crop parameters are read from *.END file
      endif

! --- set crop height and cropfactor
      if (crop%swcf.ne.3) then
        crop%common%cf = afgen (crop%fixed%cftb,(2*magrs),crop%common%dvs)
        crop%common%ch = afgen (crop%fixed%chtb,(2*magrs),crop%common%dvs)
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
                            i, state%soilwater)
        enddo
      endif

! -      n-p-k 
      if( crop%common%flCropNut) then
        call nutrinit    (nlossl,nlossr,nlosss,                         &
     &                    nuptt,rnlv,rnst,rnrt,rnso,                    &
     &                    rnldlv,rnldst,rnldrt,                         &
     &                    nfixtt,cw_nni,NLOSSLDeceasedLvToSoil)
      endif


      return

      case (2)

! === calculate potential rate and state variables =====================

! --- rates of change of the crop variables ----------------------------

! --- increase in temperature sum
      dtsum = afgen (crop%common%dtsmtb,30,atmo%Tav)

! --- phenological development rate for potential AND actual crops
      if (swsoybean.eq.0) then
! ---   standard crops
        if (crop%common%dvs.lt.1.0d0) then     
! ---     vegetative phase
          dvred = 1.0d0
          vernfac = 1.0d0
          vernrate = 0.0d0
          if (crop%wofost%idsl.ge.1) then
             dvred = max(0.0d0, min(1.0d0, (atmo%daylp - crop%wofost%dlc) / (crop%wofost%dlo - crop%wofost%dlc)))
          endif
          if (crop%wofost%idsl.eq.2) then
!            vernalisation rate,based on routines from pyWofost (Allard de Wit, 2015)
             if(.not.flvernalised) then
                if(crop%common%dvs.lt.cw_verndvs) then
                   vernrate = afgen (cw_vernrtb,30,atmo%Tav)
                   r = (vern - cw_vernbase) / (cw_vernsat - cw_vernbase)
                   vernfac = interpol(0.0d0,1.0d0,r)
                else
                   flvernalised = .true.
                endif
             endif
          endif
          dvr = vernfac * dvred*dtsum/crop%common%tsumea
        else
! ---     generative phase
          dvr = dtsum/crop%common%tsumam
        endif    

      else if (swsoybean.eq.1) then
! ---   soybean
        call mgtemprf(atmo%Tav,toptdvr,tmindvr,tmaxdvr,rfmgtemp)
        call mgphotoprf(mg,time%daynr,state%cfg%meteo%lat,popt,pcrt,flphenodayl,rfmgphotop)
        if (crop%common%dvs.lt.1.0d0) then 
! ---     vegetative phase
          if(flrfphotoveg) then
             dvr = dvrmax1 * rfmgphotop * rfmgtemp
          else
             dvr = dvrmax1 * rfmgtemp
          endif
        else
! ---     generative phase
          dvr = dvrmax2 * rfmgphotop * rfmgtemp
        endif
      endif

!     adjust development stage for realistic TSUM
      if (crop%common%dvs.ge.1.d0 .and. (.not. crop%wofost%flanthesis)) then
        crop%wofost%flanthesis = .true.
        crop%common%dvs = 1.0d0
      end if
      
! === daily dry matter production 

      gasspot = crop%wofost%pgasspot

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
      if(crop%wofost%swbulb) then
        rmrespot = (crop%common%rmr*crop%wofost%wrtpot+crop%common%rml*crop%wofost%wlvpot+crop%common%rms*crop%wofost%wstpot+crop%common%rms*crop%wofost%wblpot+ &
     &           crop%common%rmo*crop%wofost%wsopot)* afgen(crop%common%rfsetb,30,crop%common%dvs)
      else
        rmrespot = (crop%common%rmr*crop%wofost%wrtpot+crop%common%rml*crop%wofost%wlvpot+crop%common%rms*crop%wofost%wstpot+crop%common%rmo*crop%wofost%wsopot)* &
     &           afgen(crop%common%rfsetb,30,crop%common%dvs)
      endif
      teff = crop%common%q10**((atmo%Tav-25.0d0)/10.0d0)
      mrespot = dmin1(gasspot,rmrespot*teff)
      asrcpot = gasspot - mrespot

! --- partitioning factors
      fr = afgen(crop%common%frtb,30,crop%common%dvs)
      fl = afgen(crop%common%fltb,30,crop%common%dvs)
      fs = afgen(crop%common%fstb,30,crop%common%dvs)
      fo = afgen(crop%common%fotb,30,crop%common%dvs)
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
         fbl = afgen(crop%common%fbltb,30,crop%common%dvs)
      endif
! --- check on partitioning
      call chckprt(crop%common%dvs,fr,fl,fs,fo,fbl)    

! --- conversion factor 
      if(crop%wofost%swbulb) then
!       only for bulb crops (tulips etc..)
        cvf = 1.0d0/((fl/crop%common%cvl+fs/crop%common%cvs+fbl/crop%common%cvs+fo/crop%common%cvo)*(1.0d0-fr)+fr/crop%common%cvr)
      else
        cvf = 1.0d0/((fl/crop%common%cvl+fs/crop%common%cvs+fo/crop%common%cvo)*(1.0d0-fr)+fr/crop%common%cvr)
      endif
      dmipot = cvf*asrcpot
! --- check on carbon balance
      call chckcbl(crop%common%dvs,cvf,dmipot,fr,fl,fs,fbl,fo,gasspot,mrespot,      &
     &             ccheck)

! == = growth rate by plant organ

! --- growth rate roots and aerial parts
      admipot = (1.0d0-fr)*dmipot
      grrtpot = fr*dmipot
      ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
      ! growth of the roots is balanced by the death of root tissue
      if (crop%common%swrd.eq.3 .and. crop%wofost%wrtpot.gt.crop%common%wrtmax) then
        drrtpot = grrtpot
        drrtpot = max(drrtpot,crop%wofost%wrtpot*afgen (crop%common%rdrrtb,30,crop%common%dvs))
      else  
        drrtpot = crop%wofost%wrtpot*afgen (crop%common%rdrrtb,30,crop%common%dvs)
      endif  
      gwrtpot = grrtpot - drrtpot

! --- weight of new leaves
      grlvpot = fl*admipot

! --- death of leaves due to water stress or high lai
      laicr = 3.2d0/crop%kdif
      dslvpot = crop%wofost%wlvpot*max(0.0d0,min(0.03d0,0.03d0*                     &
     &           (crop%common%laipot-laicr)/laicr))

! --- death of leaves due to exceeding life span:

! --- first: leaf death due to water stress or high lai is imposed 
! ---        on array until no more leaves have to die or all leaves
! ---        are gone

      restpot = dslvpot*delt
      i1 = crop%common%ilvoldpot

      do while (restpot.gt.crop%common%lvpot(max(i1,1)).and.i1.ge.1)
        restpot = restpot - crop%common%lvpot(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalvpot = 0.0d0
      if (crop%common%lvagepot(max(i1,1)).gt.crop%common%span .and. restpot.gt.0.0d0            &
     &                   .and.i1.ge.1) then
        dalvpot = crop%common%lvpot(i1) - restpot
        restpot = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.crop%common%lvagepot(max(i1,1)).gt.crop%common%span)
        dalvpot = dalvpot+crop%common%lvpot(i1)
        i1 = i1-1
      enddo

      dalvpot = dalvpot/delt

! --- finally: calculate total death rate leaves
      drlvpot = dslvpot + dalvpot

! --- physiologic ageing of leaves per time step
      fysdel = max (0.0d0,(atmo%Tav-crop%common%tbase)/(35.0d0-crop%common%tbase))

! --- specific leaf area valid for current timestep
      slatpot = afgen (crop%common%slatb,30,crop%common%dvs)

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
      if (crop%common%laiexppot.lt.6.0d0) then
        dteff = max (0.0d0,atmo%Tav-crop%common%tbase)
! ---   increase in leaf area during exponential growth
        crop%common%glaiexpot = crop%common%laiexppot*crop%common%rgrlai*dteff
! ---   source-limited increase in leaf area
        glasolpot = grlvpot*slatpot
! ---   actual increase is determined by lowest value
        glapot = min (crop%common%glaiexpot,glasolpot)
! ---   slat will be modified in case gla equals glaiex
        if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
      endif  

! --- growth rate stems
      grstpot = fs*admipot
! --- death rate stems
      drstpot = afgen (crop%common%rdrstb,30,crop%common%dvs)*crop%wofost%wstpot
! --- net growth rate stems
      gwstpot = grstpot - drstpot

! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
! --    growth rate flowers
        grblpot = fbl*admipot
        if(crop%common%dvs.ge.1.0d0) then
          grblpot = 0.0d0
          drblpot = crop%wofost%wblpot/delt
        endif
        gwblpot = grblpot - drblpot
      endif

! --- growth rate storage organs
      gwsopot = fo*admipot

! ----integrals of the crop --------------------------------------------

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone

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

! --- leaves older than span die
      do while (crop%common%lvagepot(max(i1,1)).gt.crop%common%span.and.i1.ge.1)
        crop%common%lvpot(i1) = 0.0d0
        i1 = i1-1
      enddo

! --- oldest class with leaves
      crop%common%ilvoldpot = i1

! --- shifting of contents, updating of physiological age
      do i1 = crop%common%ilvoldpot,1,-1
        crop%common%lvpot(i1+1) = crop%common%lvpot(i1)
        crop%common%slapot(i1+1) = crop%common%slapot(i1)
        crop%common%lvagepot(i1+1) = crop%common%lvagepot(i1)+fysdel*delt
      enddo
      crop%common%ilvoldpot = crop%common%ilvoldpot + 1

! --- new leaves in class 1
      crop%common%lvpot(1) = grlvpot*delt
      crop%common%slapot(1) = slatpot
      crop%common%lvagepot(1) = 0.0d0 

! --- calculation of new leaf area and weight
      lasumpot = 0.0d0
      crop%wofost%wlvpot = 0.0d0
      do i1 = 1,crop%common%ilvoldpot
        lasumpot = lasumpot + crop%common%lvpot(i1)*crop%common%slapot(i1)
        crop%wofost%wlvpot = crop%wofost%wlvpot + crop%common%lvpot(i1)
      enddo

! --- leaf area index in case of exponential growth
      crop%common%laiexppot = crop%common%laiexppot+crop%common%glaiexpot*delt

! --- dry weight of living plant organs
      crop%wofost%wrtpot = crop%wofost%wrtpot + gwrtpot*delt
      crop%wofost%wstpot = crop%wofost%wstpot + gwstpot*delt
      crop%wofost%wsopot = crop%wofost%wsopot + gwsopot*delt
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        crop%wofost%wblpot = crop%wofost%wblpot + gwblpot*delt
      endif

! --- total above ground biomass
      tadwpot = crop%wofost%wlvpot + crop%wofost%wstpot + crop%wofost%wsopot
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
         tadwpot = tadwpot + crop%wofost%wblpot
      endif

! --- dry weight of dead plant organs (roots,leaves & stems)
      crop%wofost%dwrtpot = crop%wofost%dwrtpot + drrtpot*delt
      crop%wofost%dwlvpot = crop%wofost%dwlvpot + drlvpot*delt
      crop%wofost%dwstpot = crop%wofost%dwstpot + drstpot*delt
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        crop%wofost%dwblpot = crop%wofost%dwblpot + drblpot*delt
      endif

! --- dry weight of dead and living plant organs
      twlvpot = crop%wofost%wlvpot + crop%wofost%dwlvpot
      twstpot = crop%wofost%wstpot + crop%wofost%dwstpot
      crop%wofost%cwdmpot = twlvpot + twstpot + crop%wofost%wsopot
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        twblpot = crop%wofost%wblpot + crop%wofost%dwblpot
        crop%wofost%cwdmpot = crop%wofost%cwdmpot + twblpot
      endif

! --- total gross assimilation and maintenance respiration
      gasstpot = gasspot + gasstpot
      mrestpot = mrespot + mrestpot

! --- leaf area index
      crop%common%laipot = lasumpot + crop%common%ssa*crop%wofost%wstpot + crop%common%spa*crop%wofost%wsopot
!     prevent immediate lai reduction at emergence
!     KRO-BOO-20160403: suppressed because deviates from Wofost
!      laipot = max(laipot, laiem)

! --- vernalisation state (d)
      if(crop%wofost%idsl.eq.2) then
          vern = vern + vernrate
          if(.not.flvernalised .and. vern.ge.cw_vernsat) then
              flvernalised = .true.
          else
              if(flvernalised .and. vern.lt.cw_vernsat) then
                 write(messag,'(2a,i6,2a)') ' critical DVS,',           &
     &           ' for vernalised reached, day = ', time%daycum,          &
     &           ' but vernalisation requirements not yet fulfilled ',  &
     &           ' forcing vernalization now'
                 call log_warn('wofost', messag)
              endif
          endif
      endif


      return

      case (3)

! === calculate actual rate and state variables =====================
! === with optional calculation of (crop%common%flCropNut) water AND nutrient stress

! --- rates of change of the crop variables ----------------------------

! --- water stress reduction of pgass to gass
      if(dabs(atmo%ptra).lt.nihil) then
        crop%common%reltr = 1.0d0
      else
        crop%common%reltr = max(0.0d0,min(1.0d0,soil%tra/atmo%ptra))
      endif

! --- nitrogen stress reduction of pgass to gass
      if (crop%common%flCropNut) then
        crop%common%reltr = min(crop%common%reltr,cw_fstr)
        cw_fstr  = crop%common%reltr
      end if
      gass = crop%wofost%pgass * crop%common%reltr

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
! --  only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        rmres = (crop%common%rmr*crop%wofost%wrt+crop%common%rml*crop%wofost%wlv+crop%common%rms*crop%wofost%wst+crop%common%rms*crop%wofost%wbl+crop%common%rmo*crop%wofost%wso)* &
     &           afgen(crop%common%rfsetb,30,crop%common%dvs)
      else
        rmres = (crop%common%rmr*crop%wofost%wrt+crop%common%rml*crop%wofost%wlv+crop%common%rms*crop%wofost%wst+crop%common%rmo*crop%wofost%wso)*afgen(crop%common%rfsetb,30,crop%common%dvs)
      endif
      
      mres = dmin1(gass,rmres*teff)
      asrc = gass-mres

! --- partitioning factors
      fr = afgen(crop%common%frtb,30,crop%common%dvs)
      fl = afgen(crop%common%fltb,30,crop%common%dvs)
      fs = afgen(crop%common%fstb,30,crop%common%dvs)
      fo = afgen(crop%common%fotb,30,crop%common%dvs)
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
         fbl = afgen(crop%common%fbltb,30,crop%common%dvs)
      endif
! --- check on partitioning
      call chckprt(crop%common%dvs,fr,fl,fs,fo,fbl)    

      if( crop%common%flCropNut) then
!********************************************************************         
!         partitioning correction as influenced by water and N stress
!         Note: the partioning depends only on the Nitrogen stress,
!         not on the P and K stress. Personal communication Joost Wolf        
!         added IS
!******************************************************************** 
          CALL SUBPAR (crop%common%reltr,cw_NPART,cw_NNI,FR,FL,FS,FO)
      endif
      
! --- conversion factor 
      if(crop%wofost%swbulb) then
!       only for bulb crops (tulips etc..)
        cvf = 1.0d0/((fl/crop%common%cvl+fs/crop%common%cvs+fbl/crop%common%cvs+fo/crop%common%cvo)*(1.0d0-fr)+fr/crop%common%cvr)
      else
        cvf = 1.0d0/((fl/crop%common%cvl+fs/crop%common%cvs+fo/crop%common%cvo)*(1.0d0-fr)+fr/crop%common%cvr)
      endif
! --- dry matter increase
      dmi = cvf*asrc
! --- check on carbon balance
      call chckcbl(crop%common%dvs,cvf,dmi,fr,fl,fs,fbl,fo,gass,mres,ccheck)

! --- growth rate by plant organ

! --- growth rate roots and aerial parts
      call relgrwt(dmi,fr,fl,fs,fo,grrt,grlv,grst,grso,admi)
      if (crop%common%swrd.eq.3 .and. soil%flWrtNonox) grrt = 0.d0

! --- death of leaves due to water stress or high lai or nitrogen stress
      call deaths(crop%common%flCropNut,crop%wofost%wlv,crop%kdif,crop%lai,cw_NNI,crop%common%perdl,cw_rdrns,crop%common%reltr,dslv)

! --- death of leaves due to exceeding life span:
      call deatha(dslv,delt,crop%common%ilvold,crop%common%lv,crop%common%lvage,crop%common%span,i1,dalv)

! --- death rate leaves as result of death due to water stress or high lai and 
!                                    death due to exceeding life span
      drlv = dslv+dalv

! --- death rate stems
      drst = crop%wofost%wst * afgen (crop%common%rdrstb,30,crop%common%dvs)

! --- death rate roots
      ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
      ! growth of the roots is balanced by the death of root tissue
      if (crop%common%swrd.eq.3 .and. crop%wofost%wrt.gt.crop%common%wrtmax) then
        drrt = grrt
        drrt = max(drrt,crop%wofost%wrt*afgen (crop%common%rdrrtb,30,crop%common%dvs))
      else  
        drrt = crop%wofost%wrt*afgen (crop%common%rdrrtb,30,crop%common%dvs)
      endif  

! --- net growth rate stems, roots, storage organs
      gwst = grst - drst
      crop%wofost%gwrt = grrt - drrt
      drso = 0.0d0    ! death rate of storage organs is assumed to be 0
      gwso = grso - drso

! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
! --    growth rate flowers
        grbl = fbl*admi
        if(crop%common%dvs.ge.1.0d0) then
          grbl = 0.0d0
          drbl = crop%wofost%wbl/delt
        endif
        gwbl = grbl - drbl
      endif

! --- specific leaf area valid for current timestep
      if(crop%common%flCropNut) then
!       nutrient and water stress
        slat = afgen (crop%common%slatb,30,crop%common%dvs)*EXP(-cw_NSLA * (1.0d0-cw_NNI))
      else
        slat = afgen (crop%common%slatb,30,crop%common%dvs)
      endif
!
!     Do not allow slat higher than slatpot; slatpot can be limited by exponential growth
      slat = min(slat,slatpot)   ! pvw

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
! --- cw_FSTR is actual stress: water and nutrient 
      Fstress = crop%common%reltr
      if (crop%common%flCropNut) then
         Fstress = cw_FSTR
         if ((crop%common%dvs .LT. 0.2d0).AND.(crop%lai .LT. 0.75d0)) then
           Fstress = crop%common%reltr * EXP(-cw_NLAI* (1.0d0 - cw_NNI))
         endif
      endif
      call GLAI(Fstress,crop%common%LAIEXP,crop%common%GLAIEX,atmo%Tav,crop%common%tbase,crop%common%rgrlai,GRLV,SLAT,GLA)


! ---- UPDATE STATES: integrals of the crop --------------------------------------------

! --- phenological development stage
      crop%common%dvs = min(crop%common%dvs+dvr*delt,crop%common%dvsend)
      crop%common%tsum = crop%common%tsum + dtsum*delt

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone
      call lvdth(delt,dslv,crop%common%span,crop%common%ilvold,crop%common%lvage,crop%common%lv,i1)

! --- oldest class with leaves
      crop%common%ilvold = i1

! --- shifting of contents, updating of physiological age
      call lvshft(delt,fysdel,crop%common%ilvold,grlv,slat,crop%common%lv,crop%common%lvage,crop%common%sla)

      crop%common%ilvold = crop%common%ilvold+1

! --- calculation of new leaf area and weight
      call lvwgli(crop%common%ilvold,crop%common%lv,crop%common%sla,lasum,crop%wofost%wlv)

! --- leaf area index in case of exponential growth
      crop%common%laiexp = crop%common%laiexp+crop%common%glaiex*delt
      
! --- dry weight of living plant organs
      crop%wofost%wrt = crop%wofost%wrt+crop%wofost%gwrt*delt
      crop%wofost%wst = crop%wofost%wst+gwst*delt
      crop%wofost%wso = crop%wofost%wso+gwso*delt
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        crop%wofost%wbl = crop%wofost%wbl + gwbl*delt
      endif

! --- total above ground biomass
      tadw = crop%wofost%wlv+crop%wofost%wst+crop%wofost%wso
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        tadw = tadw + crop%wofost%wbl
      endif

! --- dry weight of dead plant organs (roots,leaves & stems)
      crop%wofost%dwrt = crop%wofost%dwrt + drrt*delt
      crop%wofost%dwlv = crop%wofost%dwlv + drlv*delt
      crop%wofost%dwst = crop%wofost%dwst + drst*delt
      crop%wofost%dwso = crop%wofost%dwso + drso*delt   ! dummy, because drso is assumed to be 0
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        crop%wofost%dwbl = crop%wofost%dwbl + drbl*delt
      endif

!     split dwlv
      idwlvCrop = (1.0d0-FraDeceasedLvToSoil) * drlv*delt
      idwlvSoil = FraDeceasedLvToSoil * drlv*delt
      crop%wofost%dwlvCrop = crop%wofost%dwlvCrop + idwlvCrop
      crop%wofost%dwlvSoil = crop%wofost%dwlvSoil + idwlvSoil
      
! --- dry weight of dead and living plant organs
!     twrt = wrt+dwrt
      twlv = crop%wofost%wlv+crop%wofost%dwlv
      twst = crop%wofost%wst+crop%wofost%dwst
      crop%wofost%cwdm = twlv+twst+crop%wofost%wso
! --- only for bulb crops (tulips etc..)
      if(crop%wofost%swbulb) then
        twbl = crop%wofost%wbl + crop%wofost%dwbl
        crop%wofost%cwdm = crop%wofost%cwdm + twbl
      endif

! --- total gross assimilation and maintenance respiration
      gasst = gass + gasst
      mrest = mres + mrest

! --- leaf area index
      crop%lai = lasum+crop%common%ssa*crop%wofost%wst+crop%common%spa*crop%wofost%wso
! --- determine maximum lai
      crop%common%laimax = max (crop%lai,crop%common%laimax)
! --- determine minimum lai to prevent dying straight after 
!       emergence when growth is slowed down due to low temperature
!     KRO-BOO-20160403: suppressed because deviates from Wofost
!      lai = max(lai, laiem)

! --- update state variables for nutrient stress

!     Calling the subroutine for N losses of leaves, roots and stem storage
!     organs (kg N ha-1 d-1)
      if(crop%common%flCropNut)then

!        Calling the subroutines for N demand of leaves, roots and stem storage
!        organs (kg N ha-1 d-1)
         CALL NDEMND(crop%wofost%wlv,crop%wofost%wst,crop%wofost%wrt,crop%wofost%wso,cw_NMAXLV,cw_NMAXST,                     &
     &                               cw_NMAXRT,cw_NMAXSO,cw_ANLV,cw_ANST,ANRT,ANSO, &
     &                               cw_tcnt,NDEML,NDEMS,NDEMR,NDEMSO)

!        Total N demand (kg N ha-1)

         NDEMTO = MAX (0.0d0,(NDEML + NDEMS + NDEMR))

!        Nutrient uptake limiting factor (-) at low moisture conditions in the
!        rooted soil layer before anthesis. After anthesis/cw_dvsnlt there is no
!        nutrient uptake from the soil
         NLIMIT = INSW(crop%common%dvs-cw_dvsnlt,INSW(crop%common%reltr-0.01d0,0.0d0,1.0d0),0.d0)
         NdemandSoil = (1.d0-cw_NFIXF) * NDEMTO * NLIMIT
         NdemandBioFix =  cw_NFIXF * NDEMTO * NLIMIT

      end if


      return

      case (4)

      if(crop%common%flCropNut) then

!        Total N uptake (kg N ha-1 d-1) from soil and by biological N fixation         
         NUPTR = (MAX(0.d0, MIN(NdemandSoil, NsupplySoil) ))/DELT
         NFIXTR = (MAX(0.d0, NdemandBioFix))/DELT

!        Calling the subroutine to estimate the translocatable nutrients in leaves, stem, roots and
!        storage organs (kg N ha-1)
         CALL NTRLOC(cw_ANLV,cw_ANST,ANRT,crop%wofost%wlv,crop%wofost%wst,crop%wofost%wrt,cw_RNFLV,cw_RNFST,cw_RNFRT,      &
     &                  cw_fntrt,ATNLV,ATNST,ATNRT,ATN)

!        N supply to the storage organs (kg N ha-1 d-1)      
         NSUPSO = INSW (crop%common%dvs-cw_dvsnt,0.0d0,ATN/cw_tcnt)

!        Rate of N uptake in grains (kg N ha-1 d-1)
         RNSO =  MIN (NDEMSO,NSUPSO)

!        Calling the subroutine to calculate nutrient translocation from leaves, stem, and roots (kg N ha-1 d-1)
         CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT)

!        Calling the subroutine to compute the partitioning of the total
!        nutrient uptake rate (NUPTR) over the leaves, stem and roots (kg N ha-1 d-1)
         CALL RNUSUB(NDEML,NDEMS,NDEMR,NUPTR,                           &
     &                  NFIXTR,NDEMTO, RNULV,RNUST,RNURT)
        
!        Calling routine to calculate nutrient losses due to dying leaves, stems           
!        and roots (kg N ha-1 d-1)    
         CALL RNLD(DRLV,DRRT,DRST,cw_RNFLV,cw_RNFRT,cw_RNFST,                    &
     &                RNLDLV,RNLDRT,RNLDST)


!   - ---Rate of change of N in crop organs   
         RNLV = RNULV - RNTLV - RNLDLV
         RNST = RNUST - RNTST - RNLDST
         RNRT = RNURT - RNTRT - RNLDRT

!-----   Total N  uptake by crop over time (kg N ha-1) from soil and by biological fixation
         NUPTT = NUPTT + NUPTR*DELT
         NFIXTT= NFIXTT+ NFIXTR*DELT

!-----   Actual N amount in various living organs and total living N amount(kg N ha-1)
         cw_ANLV =  max(0.0d0, (cw_ANLV + RNLV*DELT) )
         cw_ANST =  max(0.0d0, (cw_ANST + RNST*DELT) )
         ANRT =  max(0.0d0, (ANRT + RNRT*DELT) )
         ANSO =  ANSO + RNSO*DELT
!!!         NLIVT=  cw_ANLV + cw_ANST + ANRT + ANSO

!-----   N losses from leaves, roots and stems due to senescence and total N loss (kg N ha-1)
         NLOSSL =  NLOSSL + RNLDLV*DELT
         NLOSSR =  NLOSSR + RNLDRT*DELT
         NLOSSS =  NLOSSS + RNLDST*DELT
!!!         NLOSST =  NLOSSL + NLOSSR + NLOSSS

!----    total N  in living and dead roots
!!!         NROOT= ANRT + NLOSSR

!       increment values of dead weight of plant organs,
!       to be used in the soil nutrient submodel
         idwrt = drrt*delt
         idwlv = idwlvSoil
         idwst = 0.0d0
         idwso = 0.0d0
         iNLOSSR =  RNLDRT*DELT
         iNLOSSL =  FraDeceasedLvToSoil * RNLDLV*DELT
         NLOSSLDeceasedLvToSoil =  NLOSSLDeceasedLvToSoil + iNLOSSL
         NLOSSL =  NLOSSL - iNLOSSL
         iNLOSSS =  0.0d0
         iNLOSSO =  0.0d0
         HarLosOrm_rt = 0.0d0; HarLosOrm_lv = 0.0d0; HarLosOrm_st = 0.0d0
         HarLosOrm_dwlv = 0.0d0; HarLosOrm_dwst = 0.0d0
         HarLosOrm_so = 0.0d0; crop%common%HarLosOrm_tot = 0.0d0 
         HarLosNit_rt = 0.0d0; HarLosNit_lv = 0.0d0 
         HarLosNit_st = 0.0d0; HarLosNit_so = 0.0d0 
!!         HarLosNit_dwrt = 0.0d0; HarLosOrm_dwrt = 0.0d0
         HarLosNit_dwst = 0.0d0; HarLosNit_dwso = 0.0d0; HarLosNit_dwlv = 0.0d0 
!        during the last day of the crop period: add the weight of living roots 
!        to the dead roots and reset living weight to zero
         if (crop%common%flHarvestDay .or. (crop%common%dvs.ge.crop%common%dvsend) .or. &
     &                 dabs(time%t1900-1.0d0-crop_cfg%rotation_end(crop%common%icrop)).lt.1.0d-3 ) then
            HarLosOrm_rt = crop%wofost%wrt
            HarLosOrm_dwlv =  cw_fraharlosorm_lv * crop%wofost%dwlv
            HarLosOrm_lv   = cw_fraharlosorm_lv * crop%wofost%wlv + HarLosOrm_dwlv
            HarLosOrm_dwst =  cw_fraharlosorm_st * crop%wofost%dwst
            HarLosOrm_st   = cw_fraharlosorm_st * crop%wofost%wst + HarLosOrm_dwst
            HarLosOrm_dwso =  cw_fraharlosorm_so * crop%wofost%dwso
            HarLosOrm_so   = cw_fraharlosorm_so * crop%wofost%wso + HarLosOrm_dwso
            crop%common%HarLosOrm_tot = HarLosOrm_rt + cw_fraharlosorm_lv * crop%wofost%wlv +      &
     &             cw_fraharlosorm_st * crop%wofost%wst + cw_fraharlosorm_so * crop%wofost%wso
!ckro_sup_20170714 : suppressed because it will happen after harvest
!            wrt = wrt - HarLosOrm_rt
!            wlv = wlv - cw_fraharlosorm_lv * wlv
!            wst = wst - cw_fraharlosorm_st * wst
!            wso = wso - cw_fraharlosorm_so * wso
            idwrt = idwrt + HarLosOrm_rt
            idwlv = idwlv + HarLosOrm_lv
            idwst = idwst + HarLosOrm_st
            idwso = idwso + HarLosOrm_so
            HarLosNit_rt = ANRT
            HarLosNit_dwlv = cw_fraharlosorm_lv * NLOSSL 
            HarLosNit_lv = cw_fraharlosorm_lv * cw_ANLV + HarLosNit_dwlv 
            HarLosNit_dwst = cw_fraharlosorm_st * NLOSSS 
            HarLosNit_st = cw_fraharlosorm_st * cw_ANST + HarLosNit_dwst 
            HarLosNit_dwso = cw_fraharlosorm_so * 0.0d0 
            HarLosNit_so = cw_fraharlosorm_so * ANSO + HarLosNit_dwso 
            iNLOSSL = iNLOSSL + HarLosNit_lv
            iNLOSSS = iNLOSSS + HarLosNit_st
            iNLOSSO = iNLOSSO + HarLosNit_so
            iNLOSSR = iNLOSSR + HarLosNit_rt
            cw_ANLV = cw_ANLV - cw_fraharlosorm_lv * cw_ANLV
            cw_ANST = cw_ANST - cw_fraharlosorm_st * cw_ANST
            ANSO = ANSO - cw_fraharlosorm_so * ANSO
            ANRT = ANRT - HarLosNit_rt
        end if

! ----- CHECK and WRITE MASS BALANCE: dry matter of crop

!       output of OM balance1: from air to partitioning (kg/ha DM CH2O)
        call outbalcropOM1(2,state%cfg%general%pathwork,state%cfg%general%outfil,state%cfg%general%project,time%date,crop%common%daycrop,   &
     &       time%t,crop%common%dvs,crop%common%tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)

! -     OM balance2: storage difference(kg/ha DM CH2O)
        storagediff = (crop%wofost%wlv+crop%wofost%wst+crop%wofost%wso+crop%wofost%wrt) - (wlvt0+wstt0+wsot0+wrtt0)
        ombalan = storagediff - ( (grlv+grst+grso+grrt)*delt -          &
     &            (drlv+drst+drso+drrt)*delt )
!ckro_20171002 harvest losses happen after harvest and should not be
!    in this balance, therefore this commented out of source code.
!        storagediff = HarLosOrm_rt + HarLosOrm_lv + HarLosOrm_st +      &
!     &  HarLosOrm_so -(HarLosOrm_dwlv + HarLosOrm_dwst + HarLosOrm_dwso)
        if (dabs(ombalan) .ge. 1.0d0) then
           write(messag,'(a,f8.3)')                                     &
     &      ' Warning Wofost: OM balance2 not 0, OMBAL = ',ombalan
           call log_warn('wofost', messag)
!     &     ' OM balance2 not 0, simulation stopped OMBAL=',ombalan
!           call fatalerr ('wofost',messag)
        endif
!       output of OM balance2
        call outbalcropom2(2,state%cfg%general%pathwork,state%cfg%general%outfil,state%cfg%general%project,time%date,crop%common%daycrop,   &
     &         time%t,crop%common%dvs,crop%common%tsum,storagediff,crop%wofost%wlv,crop%wofost%wst,crop%wofost%wso,crop%wofost%wrt,delt,         &
     &         grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)

! ----- CHECK and WRITE MASS BALANCE: nitrogen of crop

        NBALAN =  NUPTT + NFIXTT + (ANLVI+ANSTI+ANRTI+ANSOI)            &
     &      - (cw_ANLV+cw_ANST+ANRT+ANSO) - (NLOSSL+NLOSSR+NLOSSS)            &
     &      - NLOSSLDeceasedLvToSoil                                    &
     &      - (HarLosNit_lv+HarLosNit_st+HarLosNit_so+HarLosNit_rt)     &
     &      +  HarLosNit_dwlv + HarLosNit_dwst + HarLosNit_dwso

!       output of N balance
        call outbalcropN(2,state%cfg%general%pathwork,state%cfg%general%outfil,state%cfg%general%project,time%date,crop%common%daycrop,     &
     &         time%t,crop%common%dvs,crop%common%tsum,NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,cw_ANLV,&
     &         cw_ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,cw_nni)
        IF (dabs(NBALAN) .GE. 1.0d-03) then
           write(messag,'(1a,i6,a,f8.3)') ' Nitrogen balance not 0,'//  &
     &     ' simulation stopped, day = ', time%daycum,' NBAL=',NBALAN
!*           call fatalerr ('wofost',messag)
           call log_warn('CropGrowth_Wofost', messag)
        endif
 
        if (crop%common%flHarvestDay .or. (crop%common%dvs.ge.crop%common%dvsend) .or. &
     &                 dabs(time%t1900-1.0d0-crop_cfg%rotation_end(crop%common%icrop)).lt.1.0d-3 ) then
          gwst  = 0.0d0
          crop%wofost%gwrt  = 0.0d0
          gwso  = 0.0d0
          grlv  = 0.0d0
          NdemandSoil = 0.0d0
          crop%common%HarLosOrm_tot = 0.0d0
        endif
      endif

! --- update normalized cumulative root density based on root extraction or stress (cumdens)
      if (crop%common%swrdc.eq.1) call update_rootdistribution(state)
      
! --- root extension
      if (crop%common%swrd.eq.1) then

        crop%common%rdpot = afgen (crop%common%rdtb,22,crop%common%dvs)
        crop%common%rdpot = min(crop%common%rdpot,crop%common%rdm)
        crop%common%rd    = crop%common%rdpot

      elseif (crop%common%swrd.eq.2) then

        rrpot = min (crop%common%rdm-crop%common%rdpot,crop%common%rri)
        if (fr.le.0.0d0 .or. crop%wofost%pgasspot.lt.1.0d0) rrpot = 0.0d0
        crop%common%rdpot = crop%common%rdpot + rrpot

        rr = min (crop%common%rdm-crop%common%rd,crop%common%rri)
        if (fr.le.0.0d0 .or. crop%wofost%pgass.lt.1.0d0 .or.                      &
     &      soil%flWrtNonox) rr = 0.0d0
        if (crop%common%swdmi2rd.eq.1 .and. crop%wofost%pgass.ge.1.0d0)              rr = rr * gass/crop%wofost%pgass
        crop%common%rd = crop%common%rd + rr

      elseif (crop%common%swrd.eq.3) then
        crop%common%rdpot = afgen (crop%common%rlwtb,22,crop%wofost%wrtpot)
        crop%common%rdpot = min(crop%common%rdpot,crop%common%rdm)
        crop%common%rd = afgen (crop%common%rlwtb,22,crop%wofost%wrt)
        crop%common%rd = min(crop%common%rd,crop%common%rdm)
      endif

! --- crop factor or crop height
      if (crop%swcf.ne.3) then
        crop%common%cf = afgen (crop%fixed%cftb,(2*magrs),crop%common%dvs)
        crop%common%ch = afgen (crop%fixed%chtb,(2*magrs),crop%common%dvs)
      else
        crop%common%cf        = afgen (crop%fixed%cftb,72,crop%lai)
        crop%fixed%cfeic = afgen (crop%fixed%cfeictb,72,crop%lai)
        crop%common%ch        = afgen(crop%fixed%chtb,72,crop%lai)
      endif
      ! cfeic write retired (was: crop%fixed%cfeic = cfeic)

! --- update canopy storage capacity
      if (crop%common%swinter.eq.3) then
        atmo%siccapact = 0.0d0   ! [GR-CROP 2026-05-25] swinter=3 stub-errored on TOML; siccaplai always 0
      endif

! --- update states of dry matter organs
      wlvt0 = crop%wofost%wlv
      wstt0 = crop%wofost%wst
      wsot0 = crop%wofost%wso
      wrtt0 = crop%wofost%wrt

      case default
         call fatalerr_collected ('Wofost', 'Illegal value for TASK')
      end select

      end associate  ! time%t1900, time%t, time%daynr, time%daycum, time%date => state%timecontrol [TC-10]
      return
      end subroutine wofost

! ----------------------------------------------------------------------
      subroutine outbalcropN(task,pathwork,outfil,project,date,daycrop, &
     &             t,dvs,tsum,NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV,&
     &               ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,NNI) 
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output N balance files
! ----------------------------------------------------------------------
      use error_mod, only: fatalerr_collected
      use file_io_mod, only: file_open
      implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum   !,crop%common%laipot,lai,crop%common%cf,rdpot,rd,crop%common%ch,crt0,crt1
!      real(8) cwdmpot,cwdm,wsopot,wso,wstpot,wst,wlvpot,wlv,wrtpot,wrt
      real(8) NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV
      real(8) ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,NNI
! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   nba
      character(len=160) filnam,filtext

      save    nba

      comma = ',' 
    
      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.nba'
      call file_open(nba, filnam, 'replace', 'write')
      filtext = 'output of N-balance of detailed crop growth model'
      call writehead (nba,1,filnam,filtext,project)

      write (nba,100)
 100    format ('*',/,                                                  &
     & '*             day     day      -    grC   kg/ha   kg/ha   kg/', &
     & 'ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   ', &
     & 'kg/ha   kg/ha   kg/ha  kg/ha',/,                                &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,  NUPTT, NFIXTT,  ANLVI'&
     & ,',  ANSTI,  ANRTI,  ANSOI,   ANLV,   ANST,   ANRT',             &
     & ',   ANSO, NLOSSL, NLOSSR, NLOSSS, NBALAN, NNI')
      return

      case (2)
! --- write dynamic data ----------------------------------------------------
      write (nba,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,NUPTT,comma,NFIXTT,comma,ANLVI,comma,ANSTI,comma,     &
     & ANRTI,comma,ANSOI,comma,ANLV,comma,ANST,comma,ANRT,comma,        &
     & ANSO,comma,NLOSSL,comma,NLOSSR,comma,NLOSSS,comma,NBALAN,        &
     & comma,NNI
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 15(a1,f7.2) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (nba)

      case default
         call fatalerr_collected ('OutbalCropN', 'Illegal value for TASK')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outbalcropOM1(task,pathwork,outfil,project,date,       &
     &      daycrop,t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output OM balance files 
! ----------------------------------------------------------------------
      
     use swap_constants, only: tiny, nihil
     use error_mod, only: fatalerr_collected
     use file_io_mod, only: file_open
     implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck
! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   om1
      character(len=160) filnam,filtext
      real(8)   OMroot,OMleaves,OMstems,OMstorage  !,Cccheck

      save    om1

      comma = ',' 

      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.om1'
      call file_open(om1, filnam, 'replace', 'write')
      filtext = 'output of OM1-balance (kg/ha DM increase per time'//   &
     & 'step) of detailed crop growth model'
      call writehead (om1,1,filnam,filtext,project)

      write (om1,100)
 100    format ('*',/,                                                  &
     & '*             day     day      -   grCd   kg/ha   kg/ha',       &
     & '   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha       -   kg/ha',/,    &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,   gass,   mres,OMroots'&
     & ,',OMleaves,OMstems,OMstorage, dmi,    cvf,OMcheck')
      return

      case (2)

! --- write dynamic data ----------------------------------------------------
      if(cvf.lt.nihil) then
        OMroot    = 0.0d0 
        OMleaves  = 0.0d0
        OMstems   = 0.0d0
        OMstorage = 0.0d0
      else
        OMroot    = fr*dmi/cvf
        OMleaves  = fl*(1.0d0-fr)*dmi/cvf
        OMstems   = fs*(1.0d0-fr)*dmi/cvf
        OMstorage = fo*(1.0d0-fr)*dmi/cvf
      endif  

      write (om1,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,gass,comma,mres,comma,OMroot,comma,OMleaves,comma,    &
     & OMstems,comma,OMstorage,comma,dmi,comma,cvf,comma,ccheck
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 9(a1,f7.2) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (om1)

      case default
         call fatalerr_collected ('OutbalCropOM1', 'Illegal value for TASK')
      end select

      return
      end 

      subroutine outbalcropom2(task,pathwork,outfil,project,date,       &
     &         daycrop,t,dvs,tsum,storagediff,wlv,wst,wso,wrt,      &
     &         delt,grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output OM balance files
! ----------------------------------------------------------------------
      use error_mod, only: fatalerr_collected
      use file_io_mod, only: file_open
      implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt
      real(8) grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan

! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   om2
      character(len=160) filnam,filtext

      save    om2

      comma = ',' 

      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.om2'
      call file_open(om2, filnam, 'replace', 'write')
      filtext = 'output of OM2-balance (kg/ha DM, cumulative and '//    &
     & 'increments) of detailed crop growth model'
      call writehead (om2,1,filnam,filtext,project)

      write (om2,100)
 100    format ('*',/,                                                  &
     & '*             day     day      - degday   kg/ha   kg/ha ',      &
     & '  kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha ',      &
     & '  kg/ha kg/ha/d kg/ha/d kg/ha/d kg/ha/d',/,                     &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,storagediff,wlv,',      &
     & '    wst,    wso,    wrt,   grlv,   grst,   grso,   grrt,',      &
     & '   drlv,   drst,   drso,   drrt,ombalan')

      return

      case (2)

! --- write dynamic data ----------------------------------------------------
      write (om2,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,storagediff,comma,wlv,comma,wst,comma,wso,comma,      &
     & wrt,comma,grlv*delt,comma,grst*delt,comma,grso*delt,comma,       &
     & grrt*delt,comma,drlv*delt,comma,drst*delt,comma,drso*delt,       &
     & comma,drrt*delt,comma,ombalan
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 14(a1,f7.1) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (om2)

      case default
         call fatalerr_collected ('OutbalCropOM2', 'Illegal value for TASK')
      end select

      return
      end 


      SUBROUTINE CHCKPRT(DVS,FR,FL,FS,FO,FBL)       
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Checks the partitioning factors, and interrupt in case of error
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    R8  FR        Fraction of total dry matter partitioned to the roots (-)
!       I    R8  FL        Fraction of total dry matter partitioned to the leaves (-)
!       I    R8  FS        Fraction of total dry matter partitioned to the stems (-)
!       I    R8  FO        Fraction of total dry matter partitioned to the storage organs (-)
!       I    R8  FBL       Fraction of total dry matter partitioned to the bulbs (-)
! ----------------------------------------------------------------------
      use error_mod, only: fatalerr_collected
      implicit none
! --- global
      real(8)   DVS,FR,FL,FS,FO,FBL
! --- local
      character(len=300) messag
      real(8)   FCHECK

!*     check on partitioning
      fcheck = fr+(fbl+fl+fs+fo)*(1.0d0-fr) - 1.0d0
      if (dabs (fcheck).gt.0.0001d0) then
!        write (messag,'(a,f5.2,/,3(a,g12.5),/,2(a,g12.5))')             &
        write (messag,'(a,f5.2,3(a,g12.5),2(a,g12.5))')                 &
     &      ' error in partitioning functions, dvs= ',dvs,              &
     &      ' fcheck = ',fcheck,' fr = ',fr,' fl = ',fl,                &
     &      ' fs = ',fs,' fo = ',fo
        call fatalerr_collected ('wofost',messag)
      end if
      
      return
      end
 
      subroutine chckcbl(dvs,cvf,dmi,fr,fl,fs,fbl,fo,gass,mres,ccheck)
! ----------------------------------------------------------------------
!     Last modified      : July 2017
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Checks the carbon balance, and interrupt in case of error
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    R8  DVS       development stage (-)
!       I    R8  CVF       conversion factor (-)
!       I    R8  DMI       dry matter increase (-)
!       I    R8  FBL       Fraction of total dry matter partitioned to the bulbs (-)
!       I    R8  FR        Fraction of total dry matter partitioned to the roots (-)
!       I    R8  FL        Fraction of total dry matter partitioned to the leaves (-)
!       I    R8  FS        Fraction of total dry matter partitioned to the stems (-)
!       I    R8  FO        Fraction of total dry matter partitioned to the storage organs (-)
!       I    R8  GASS      Gross assimilation (-)
!       I    R8  MRES      maintenance respiration (-)
! ----------------------------------------------------------------------
      use error_mod, only: fatalerr_collected
      implicit none
! --- global
      real(8)   dvs,cvf,dmi,fbl,fr,fl,fs,fo,gass,mres,ccheck
! --- local
      character(len=300) messag

!     check on c-balance
      ccheck = (gass-mres-(fr+(fbl+fl+fs+fo)*(1.0d0-fr))*dmi/cvf)       &
     &       /max (0.0001d0,gass)
      if (dabs (ccheck).gt.0.0001d0) then
        write (messag,'(a,i3,/,3(a,g12.5),/,a,4g12.5,/,2(a,g12.5))')    &
     &     ' carbon flows nog balanced on day ',dvs,                    &
     &     ' ccheck = ',ccheck,' gass = ',gass,' mres = ',mres,         &
     &    ' fr,fbl,l,s,o = ',fr,fbl,fl,fs,fo,' dmi = ',dmi,' dvf = ',cvf
        call fatalerr_collected ('wofost',messag)
      end if    
      return
      end

      subroutine relgrwt(dmi,fr,fl,fs,fo,grrt,grlv,grst,grso,admi)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : To calculate relative growth rate of roots, stems leaves
!                          and storage organs  
! Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DSLV     Death rate leaves                              kg ha-1 d-1 DM
!   I    R8  DMI      total dry matter increase                      kg ha-1 d-1 DM
!   I    R8  FR       Fraction of total dry matter partitioned to the roots (-)
!   I    R8  FL       Fraction of total dry matter partitioned to the leaves (-)
!   I    R8  FS       Fraction of total dry matter partitioned to the stems (-)
!   I    R8  FO       Fraction of total dry matter partitioned to the storage organs (-)
!   O    R8  GRRT     Growth rate roots                              kg ha-1 d-1 DM
!   O    R8  GRLV     Growth rate leaves                             kg ha-1 d-1 DM
!   O    R8  GRST     Growth rate stems                              kg ha-1 d-1 DM
!   O    R8  GRSO     Growth rate storage organs                     kg ha-1 d-1 DM
!   O    R8  ADMI     total above ground dry matter increase         kg ha-1 d-1 DM
! ===== ==== =======  =============================================  ==============
!   O    R8  ADMI      above ground dry matter increase (-)
      implicit none
! --- global
      real(8)   DMI,FR,FL,FS,FO,GRRT,GRLV,GRST,GRSO,ADMI
! --- local
!     save
            
      ADMI = (1.0d0-FR)*DMI
      GRRT = FR*DMI
      GRLV = FL*ADMI
      GRST = FS*ADMI
      GRSO = FO*ADMI
      
      RETURN
      END

      subroutine deaths(flcropnut,wlv,kdif,lai,NNI,perdl,rdrns,         &
     &                  reltr,dslv)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Compute the relative death rate leaves due
!                          to stress (kg DM ha-1 d-1)   
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    L   flCropNut Flag indicating simulation of nutrient stress -
!   I    R8  WLV       Dry weight of living leaves                   kg ha-1 DM
!   O    R8  reltr
!   I    R8  KDIF      Extinction coeff. for diffuse visible light   -
!   I    R8  LAI       Leaf Area Index
!   I    R8  NNI  
!   I    R8  PERDL     Max.rel.death rate of leaves due to water strs -
!   I    R8  RDRNS     Max.rel.death rate of leaves due to nitrogen strs -
!   O    R8  DSLV      Death rate leaves                             kg ha-1 d-1 DM
! ===== ==== =======  =============================================  ==============
!   -    R8  DSLV1     Death rate leaves due to water stress         kg ha-1 d-1 DM
!   -    R8  DSLV1     Death rate leaves due to self-shading         kg ha-1 d-1 DM
!   -    R8  LAICR     
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      logical   flCropNut
      real(8)   DSLV,KDIF,LAI,NNI,PERDL,RDRNS,reltr,WLV
! --- local
      real(8)   DSLV1,DSLV2,LAICR
!      save     
      
!     death rate of leaves due to water stress
      DSLV1 = WLV*(1.d0-reltr)*PERDL
      
!     death rate of leaves due high LAI
      LAICR = 3.2d0/KDIF
      DSLV2 = WLV*max(0.0d0, min(0.03d0, (0.03d0*(LAI-LAICR)/LAICR)))
      DSLV  = MAX (DSLV1, DSLV2)
      
!     death rate increase due to nutrient shortage
      IF(flCropNut .AND. NNI.LT.1.0d0) THEN
         DSLV = DSLV + WLV*RDRNS * (1.0d0-NNI)
      END IF 
      
      RETURN
      END
      
      subroutine deatha(dslv,delt,ilvold,lv,lvage,span,i1,dalv)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : To compute the relative death rate leaves due              * 
!                          to ageing (kg DM ha-1 d-1)   

!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT      
!   I    R8  REST   
!   I    I   i1vold      
!   I    R8  SPAN  
!   I    R8  LVAGE
!   I/O  I   11      
!   O    R8  DALV      Death rate leaves due to ageing (kg/ha/d DM)
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   i1,ilvold
      real(8)   DSLV,DELT,SPAN,DALV,lv(366),lvage(366)
! --- local
      real(8)   REST
            
! --- first: leaf death due to water stress or high lai is imposed on array
! ---        until no more leaves have to die or all leaves are gone

      rest = dslv*delt
      i1 = ilvold

      do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
        rest = rest-lv(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalv = 0.0d0
      if (lvage(max(i1,1)).gt.span .and. rest.gt.0.0d0 .and.i1.ge.1)then
        dalv = lv(i1)-rest
        rest = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
        dalv = dalv+lv(i1)
        i1 = i1-1
      enddo

      dalv = dalv/delt
      
      RETURN
      END      

      SUBROUTINE GLAI(Fstress,LAIEXP,GLAIEX,TEMP,TBASE,RGRLAI,GRLV,SLAT,GLA)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : exponential, sink limited leave increase
!                          and adjustment of specific leaf area of youngest 
!                          leaf class. Adjust for water and nitrogen stress
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  Fstress
!   I/O  R8  LAIEXP
!   I/O  R8  GLAIEX
!   O    R8  GLA       
!   I/O  R8  SLAT     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      real(8)   Fstress,LAIEXP,GLAIEX,TEMP,TBASE,RGRLAI,GRLV,SLAT,GLA
! --- local
      real(8)   DTEFF,GLASOL
!      save     

      IF (LAIEXP.LT.6.0d0) THEN
         DTEFF  = MAX (0.d0,TEMP-TBASE)
         GLAIEX = Fstress * LAIEXP*RGRLAI*DTEFF
!*        source-limited increase in leaf area
         GLASOL = GRLV*SLAT
!*        sink-limited increase in leaf area
         GLA    = MIN (GLAIEX, GLASOL)
!*        adjustment of specific leaf area of youngest leaf class
         IF (GRLV.GT.0.d0) SLAT = GLA/GRLV
      END IF
      
      RETURN
      END
      SUBROUTINE LVDTH(DELT,DSLV,SPAN,ILVOLD,LVAGE,LV,I1)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Impose leave death on LV array 
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT
!   I    R8  DSLV       
!   I    I   SPAN
!   I    I   ILVOLD
!   I    R8  LVAGE    
!   I/O  R8  LV     
!   O    I   I1
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold,i1
      real(8)   delt,dslv,lv(366),lvage(366),span
! --- local
      real(8)   dslvt

      save     

! --- remaining leaves
      dslvt = dslv*delt
      i1 = ilvold
      do while ((dslvt.gt.0.0d0) .and. (i1.ge.1))
        if (dslvt.ge.lv(i1)) then
          dslvt = dslvt-lv(i1)
          lv(i1) = 0.0d0
          i1 = i1-1
        else
          lv(i1) = lv(i1)-dslvt
          dslvt = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvage(max(i1,1)).gt.span.and.i1.ge.1)
        lv(i1) = 0.0d0
        i1 = i1-1
      enddo

      RETURN
      END

      SUBROUTINE LVSHFT(DELT,FYSDEL,ILVOLD,grlv,slat,LV,LVAGE,SLA)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Shift contents of LV, LVAGE, SLA tables with one day  
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT
!   I    R8  FYSDEL       
!   I    I   ILVOLD
!   I/O  R8  LV     
!   I/O  R8  LVAGE     
!   I/O  R8  SLA     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold
      real(8)   delt,fysdel,grlv,slat,lv(366),lvage(366),sla(366)
! --- local
      integer   i1

!      save     

!     Shift contents of LV, LVAGE, SLA tables with one day
      do i1 = ilvold,1,-1
        lv(i1+1) = lv(i1)
        sla(i1+1) = sla(i1)
        lvage(i1+1) = lvage(i1)+fysdel*delt
      enddo

!     new leaves in class 1
      lv(1) = grlv*delt
      sla(1) = slat
      lvage(1) = 0.0d0 
     
      return
      end
      SUBROUTINE LVWGLI(ILVOLD,LV,SLA,LASUM,WLV)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Calculate new wlv and lai 
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    I   ILVOLD
!   I    R8  LV     
!   I    R8  SLA     
!   O    R8  LASUM     
!   O    R8  WLV     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold
      real(8)   lv(366),sla(366),lasum,wlv
! --- local
      integer    i1
!      save     

! --- calculation of new leaf area and weight
      lasum = 0.0d0
      wlv = 0.0d0
      do i1 = 1,ilvold
        lasum = lasum+lv(i1)*sla(i1)
        wlv = wlv+lv(i1)
      enddo
  
      return
      end

! --- for soybean (swsoybean=1): temperature and photoperiodicity
      subroutine mgtemprf(tav,toptdvr,tmindvr,tmaxdvr,rfmgtemp)
! ----------------------------------------------------------------------
!     Last modified      : Sept 2015
!       based on routines needed from pyWofost for soybean (Allard de Wit, 2015)
!
!     Purpose            : temperature reduction factor for soybean (short day)
!       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
!       http://digitalcommons.unl.edu/agronomyfacpub/112
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  tav        
!   I    R8  toptdvr       
!   I    R8  tmindvr
!   I    R8  tmaxdvr     
!   O    R8  rfmgtemp    
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      real(8)   tav, toptdvr, tmindvr, tmaxdvr, rfmgtemp
! --- local
      real(8)   alpha, p1, p2, p3, p4

      alpha = log(2.0d0)/(log((tmaxdvr-tmindvr)/(toptdvr-tmindvr)))
      if(tav.lt.tmindvr .or. tav.gt.tmaxdvr) then
        rfmgtemp = 0.0d0
      else
        p1 = 2.0d0 * (tav - tmindvr)**alpha
        p2 = (toptdvr - tmindvr)**alpha
        p3 = (tav - tmindvr)**(2.0d0*alpha)
        p4 = (toptdvr - tmindvr)**(2.0d0*alpha)
        rfmgtemp = (p1 * p2 - p3) / p4

      endif

      return
      end

      subroutine mgphotoprf(mg,iday,lat,popt,pcrt,flphenodayl,rfmgphotop)
! ----------------------------------------------------------------------
!     Last modified      : Sept 2015
!       based on routines needed from pyWofost for soybean (Allard de Wit, 2015)
!
!     Purpose            : Photoperiod reduction factor for soybean (short day)
!       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
!       http://digitalcommons.unl.edu/agronomyfacpub/112
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  mg       maturity group 
!   I    I   iday     daynr      
!   I    R8  lat      lattitude
!   I    L   flphenodayl Flag to allow input of POPT and PCRT or using 
!                        empirical relation from Setiyono et al
!   I    R8  popt     optimal daylength for phenological developm.   hr
!   I    R8  pcrt     critical daylength for phenological developm.  hr
!   O    R8  rfmgphotop  reduction factor for photoperiodicity       -  
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   iday
      logical   flphenodayl
      real(8)   lat, popt, pcrt
      real(8)   mg, rfmgphotop
! --- local
      real(8)   alpha,m,p0,p1,p2
      real(8)   dec,daylp,pi,rad,sinld,cosld,aob
      parameter (pi=3.1415926d0, rad=0.0174533d0)

! astronomic daylength according to solar elevation angle of -0.833 day
!*     Declination and solar constant for this day
      dec = -asin(dsin(23.45d0*rad)*dcos(2.d0*pi*dble(iday+10)/365.0d0))
      sinld = dsin(rad*lat)*dsin(dec)
      cosld = dcos(rad*lat)*dcos(dec)
      aob = sinld/cosld
      if (dabs(aob).le.1.0d0) then 
         daylp = 12.0d0*(1.d0+2.d0*asin(aob)/pi)
      else
         if (aob.gt.1.0d0)  daylp = 24.0d0
         if (aob.lt.-1.0d0) daylp =  0.0d0
      endif

! First determine Popt and Pcrt based on maturity group rating
      m = 3.0d0
      if(flphenodayl) then
        continue            ! popt and pcrt  are input
      else
        popt = 12.759d0 - 0.388d0*mg - 0.058d0*mg**2
        pcrt = 27.275d0 - 0.493d0*mg - 0.066d0*mg**2
      endif
      alpha = log(2.0d0)/log(((pcrt - popt)/m) + 1.0d0)
      p0 = (pcrt - popt)/m

      if(daylp.lt.popt) then
        rfmgphotop = 1.0d0
      else
        if (daylp.gt.pcrt) then
          rfmgphotop = 0.0d0
        else
          p1 = (daylp - popt)/m + 1.0d0
          p2 = (pcrt - daylp)/(pcrt - popt)
          rfmgphotop = (p1*(p2**p0))**alpha
        endif
      endif

!      write(99,*)  iday, daylp, rfmgphotop
      
      return
      end


      end module cropwofost_runtime_mod
