! File VersionID:
!   $Id: initialize.f90 374 2018-03-21 13:12:23Z heine003 $
! ----------------------------------------------------------------------
      Subroutine Initialize
! --- Initialize all variables (except crop) in subroutine Variables

      use variables
      implicit none 

! --- time & control variables
      ! [SS-BMI2 Task 5] retired 2026-05-13 — all 15 TimeControl config fields moved to
      ! state%timecontrol, populated by config_to_variables before timecontrol_init runs (ADR 0041)
      ! flprintdt, nprintday, period, swheader, swodat, swres, swscre, dtmax, dtmin,
      ! tend, tstart, MaxIt, MaxIterTime, msteps, flMaxIterTime — deleted.
      ! [SS-TC] retired 2026-05-12 — TimeControl runtime fields default to zero in state%timecontrol (ADR 0041)
      logf               = 0
      ! [SS-GR-FINAL D1] ex_tlast zero-fill dropped — declaration retired
      outdat             = 0.0d0
      outdatint          = 0.0d0
      flSwapShared       = .false.

! --- crop variables (not crop specific, such as crop calendar)
      crp                = 0 
      icrop              = 0 
      croptype           = 0 
      cropend            = 0.0d0 
      cropstart          = 0.0d0 
      rdmax              = 0.0d0 
      flCropOutput       = .false. 
      flCropCalendar     = .false.
      flCropOpenFile     = .false.
      flCropReadFile     = .false.
      flCropEmergence    = .false.
      flCropHarvest      = .false.
      
! --- meteo variables
      ! [SS-TC] daymeteo/rainrec/swmeteo/wrecord/yearmeteo retired to state%timecontrol (ADR 0041)
      idregr             = 0
      iharvest           = 1
      ilvold             = 0
      ilvoldpot          = 0
      iseqgm             = 0
      iseqgmpot          = 0
      nmetdetail         = 0
      swdivide           = 0
      swetr              = 0
      swetsine           = 0
      swinter            = 0
      swmetdetail        = 0
      swrain             = 0
      ! [SS-ATM A-2.6] aintcdt/atmdem/caintc/cevap/cgrai/cnrai/cpeva/cptra retired to state%atmosphere
      alt                = 0.0d0
      altw               = 0.0d0
      daylp              = 0.0d0

! --- only for bulb crops (tulips etc..)
      dwbl             = 0.0d0
      dwblpot          = 0.0d0
      dwso             = 0.0d0
      
      dwlv               = 0.0d0 
      dwlvCrop           = 0.0d0 
      dwlvSoil           = 0.0d0 
      dwlvpot            = 0.0d0 
      dwrt               = 0.0d0 
      dwrtpot            = 0.0d0 
      dwst               = 0.0d0 
      dwstpot            = 0.0d0 
      ! [SS-ATM A-2.6] empreva/fprecnosnow/grai/graidt/ievap/inrai/ipeva/iptra retired to state%atmosphere
      epot               = 0.0d0
      cfevappond         = 0.0d0
      finterception      = 0.0d0
      grain              = 0.0d0
      laiexp             = 0.0d0 
      laiexppot          = 0.0d0 
      laimax             = 0.0d0 
      lat                = 0.0d0 
      ! [SS-TC] metperiod retired to state%timecontrol%metperiod (ADR 0041)
      ! [SS-ATM A-2.6] nraida/nraidt/peva/pevaday/ptra/ptraday retired to state%atmosphere
      nrain              = 0.0d0
      rad                = 0.0d0
      ! [GR-CROP C12] rainamount/rainfluxarray/raintimearray retired from variables — state%atmosphere%X
      raintab            = 0.0d0
      ! [GR-ATM C8] rh retired — state%atmosphere%rh written by meteoday
      tadw               = 0.0d0
      tadwpot            = 0.0d0
      tav                = 0.0d0
      ! [GR-ATM C8] tavd retired — state%atmosphere%tavd written by meteoday
      tmn                = 0.0d0 
      tmnr               = 0.0d0 
      tmx                = 0.0d0 
      tpot               = 0.0d0 
      ! [SS-SWC] tra retired — soilwater_init handles init via state%soilwater%tra
      ! tra                = 0.0d0
      ! [SS-TC] fletsine/flmeteodt/flmetdetail/flrainintens retired to state%timecontrol (ADR 0041)

! --- irrigation variables
      ! [SS-GR-FINAL D1] irg zero-fill dropped — declaration retired
      irrigevent         = 0
      irtype             = 0
      isua               = 0
      isuas              = 0
      nirri              = 0
      ! [SS-GR-FINAL D1] phormc zero-fill dropped — declaration retired 
      schedule           = 0 
      swirfix            = 0 
      cirrs              = 0.0d0 
      dcrit              = 0.0d0 
      ditab              = 0.0d0 
      dwatab             = 0.0d0 
      fidtab             = 0.0d0 
      gird               = 0.0d0 
      hcritab            = 0.0d0 
      ! [SS-SWC] igird/inird retired — soilwater_init handles init via state%soilwater
      ! igird              = 0.0d0
      ! inird              = 0.0d0
      irconc             = 0.0d0 
      irdate             = 0.0d0 
      irdepth            = 0.0d0 
      nird               = 0.0d0 
      raithreshold       = 0.0d0 
      rawtab             = 0.0d0 
      tawtab             = 0.0d0 
      tcritab            = 0.0d0 
      tstairrig          = 0.0d0 
      tendirrig          = 0.0d0 
      treltab            = 0.0d0 
      ! [SS-TC] flheadirg/flirrigate/flIrg1Start retired to state%timecontrol (ADR 0041)
      ! [SS-GR-FINAL D1] FlIrrigationOutput zero-fill dropped — declaration retired 

! --- soilwater variables
      MaxBackTr          = 0
      ! [SS-BMI2 Task 5] MaxIt retired 2026-05-13 — moved to state%timecontrol%MaxIt
      Itnumb             = 0
      CritDevh1Cp        = 0.0d0
      CritDevh2Cp        = 0.0d0
      CritDevPondDt      = 0.0d0
      fldumpconvcrit     = .false.
      flwarn_hc          = .true.     ! Initialize headcalc warning flag (previously SAVE variable)
      iwarn_hc           = 0          ! Initialize headcalc warning counter (previously SAVE variable)
      ! [SS-GR-FINAL D1] nstep_hc zero-fill dropped — declaration retired 
      dev_cmb            = 0          ! Initialize mass balance deviation file unit (previously SAVE in checkmassbal)
      ! [GR-BH Task 36] flksatexm retired — seeded via state%soilwater%flksatexm in soilwater_init
      ! flksatexm          = .false.
      ! [SS-SWC] fluseksatexm retired — soilwater_init handles init via state%soilwater%fluseksatexm
      ! fluseksatexm       = .false.
      ! [SS-BND B-2.7] FlRunoff retired — soilwater_init handles init via state%soilwater
      ! FlRunoff           = .false.
      ! [GR-BH Task 36] H0max/k1max/q0 retired — seeded via state%soilwater in swap_mod.f90
      ! h0max              = 0.0d0
      ! k1max              = 0.0d0
      ! q0                 = 0.0d0
      ! [SS-GR-FINAL D1] afo/aun/bal/blc/bma zero-fills dropped — declarations retired 
      ! [SS-SWC] bpegwl retired — soilwater_init handles init via state%soilwater%bpegwl
      ! bpegwl             = 0
      botcom             = 0 
      dra                = 0 
      dramet             = 0 
      inc                = 0 
      ! [SS-SWC] indeks retired — soilwater_init handles init via state%soilwater%indeks
      ! indeks             = 0
      ipos               = 0 
      isoillay           = 0 
      ! [GR-BH Task 35] layer retired to state%mesh%layer
      ! [SS-BMI2 Task 5] msteps retired 2026-05-13 — moved to state%timecontrol%msteps
      ncomp              = 0 
      nhead              = 0 
      ! [SS-SWC] nodgwl retired — soilwater_init handles init via state%soilwater%nodgwl
      ! nodgwl             = 0
      nod1lay            = 0 
      ! [SS-HEAT] Task 9: nodfrostbot retired to state%heat%nodfrostbot
      ! nodfrostbot        = 0
      ! [SS-SWC] npegwl retired — soilwater_init handles init via state%soilwater%npegwl
      ! npegwl             = 0
      ! [GR-BH Task 37] nrlevs retired to state%drainage%nrlevs
      nsublay            = 0
      numbit             = 0 
      numlay             = 0 
      ! [GR-BH Task 35] numnod retired to state%mesh%numnod
      numnodnew          = 0 
      numtab             = 0 
      numtablay          = 0 
      ! [SS-GR-FINAL D1] str zero-fill dropped — declaration retired
      sw2                = 0
      sw3                = 0
      sw4                = 0
      ! [SS-GR-FINAL D1] swafo/swaun/swini/swbal/swblc/swwba/swsba zero-fills dropped — declarations retired
      ! [GR-BH Task 36] swbotb retired — seeded via state%soilwater%swbotb_runtime in swap_mod.f90
      ! swbotb             = 0
      swbotbc            = 0 
      swbotb3Impl        = 0 
      SwBotb3ResVert     = 0 
      swcfbs             = 0 
      swdiscrvert        = 0 
      ! [GR-BH Task 37] swdivd retired to state%drainage%swdivd
      ! [GR-BH Task 37] swdivdinf retired to state%drainage%swdivdinf
      swdislay           = 0
      swtopdislay        = 0
      swdra              = 0 
      swfrost            = 0 
      swhyst             = 0 
      swinco             = 0 
      swkmean            = 0 
      swkimpl            = 0 
      ! [SS-GR-FINAL D5] swoutputmodflow zero-fill dropped — declaration retired
      swpondmx           = 0 
      swqhbot            = 0 
      swredu             = 0 
      swsophy            = 0
      ! [SS-GR-CROPRT C1] swend zero-fill dropped — global + state field retired (ADR 0009: always 0)
      ! [SS-GR-FINAL D1] swvap/vap/wba zero-fills dropped — declarations retired
      aqamp              = 0.0d0 
      aqave              = 0.0d0 
      aqper              = 0.0d0 
      aqtmax             = 0.0d0 
      basegw             = 0.0d0 
      bdens              = 0.0d0 
      cfbs               = 0.0d0 
      ! [SS-SWC] cgird/cinund/cnird retired — soilwater_init handles init via state%soilwater
      ! cgird              = 0.0d0
      ! cinund             = 0.0d0
      ! cnird              = 0.0d0
      ! [GR-BH Task 36] cofani retired — seeded via state%soilwater%cofani in swap_mod.f90
      ! cofani             = 0.0d0
      ! [SS-SWC] cofgen retired — soilwater_init handles init via state%soilwater%cofgen
      ! cofgen             = 0.0d0
      cofqha             = 0.0d0 
      cofqhb             = 0.0d0 
      cofred             = 0.0d0 
      ! [SS-SWC] cqbot/cqbotdo/cqbotup retired — soilwater_init handles init via state%soilwater
      ! cqbot              = 0.0d0
      ! cqbotdo            = 0.0d0
      ! cqbotup            = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: cqdra/cqdrain/cqdrainin/cqdrainout removed (now in state%surfacewater).
      ! cqdra=0, cqdrain=0, cqdrainin=0, cqdrainout=0 — initialized in surfacewater_state_t defaults.
      ! [SS-SWC] cqprai/cqssdi/cqtdo/cqtup retired — soilwater_init handles init via state%soilwater
      ! cqprai             = 0.0d0
      ! cqssdi             = 0.0d0
      ! cqtdo              = 0.0d0
      ! cqtup              = 0.0d0
      CritDevMasBal      = 0.0d0 
      ! [SS-SWC] crunoff/crunoffCN/crunon retired — soilwater_init handles init via state%soilwater
      ! crunoff            = 0.0d0
      ! crunoffCN          = 0.0d0
      ! crunon             = 0.0d0
      cseeptab           = 0.0d0
      ! [SS-BND B-2.7] deepgw retired — soilwater_init handles init via state%soilwater
      ! deepgw             = 0.0d0
      ! [SS-SWC] dimoca retired — soilwater_init handles init via state%soilwater%dimoca
      ! dimoca             = 0.0d0
      ! [GR-BH Task 35] disnod retired to state%mesh%disnod
      ! drainl             = 0.0d0  ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      drares             = 0.0d0 
      ! [GR-BH Task 35] dz retired to state%mesh%dz
      dznew              = 0.0d0 
      entres             = 0.0d0 
      es0                = 0.0d0 
      et0                = 0.0d0 
      ew0                = 0.0d0 
      ! [SS-SWC] evp retired — soilwater_init handles init via state%soilwater%evp
      ! evp                = 0.0d0
      ! [GR-BH Task 37] FacDpthInf retired to state%drainage%FacDpthInf
      ftopdislay         = 0.0d0
      geofac             = 0.0d0 
      ! [SS-SWC] gwl retired — soilwater_init handles init via state%soilwater%gwl
      ! gwl                = 0.0d0
      gwlconv            = 0.0d0 
      gwli               = 0.0d0 
      ! [SS-BND B-2.7] gwlinp retired — soilwater_init handles init via state%soilwater
      ! gwlinp             = 0.0d0
      ! [SS-SWC] gwlm1 retired — soilwater_init handles init via state%soilwater%gwlm1
      ! gwlm1              = 0.0d0
      gwltab             = 0.0d0 
      ! [SS-SWC] h retired — soilwater_init handles init via state%soilwater%h
      ! h                  = 0.0d0
      h_enpr             = 0.0d0 
      haqtab             = 0.0d0 
      ! [SS-SWC] hatm retired — soilwater_init handles init via state%soilwater%hatm
      ! hatm               = 0.0d0
      ! [SS-BND B-2.7] hbot retired — soilwater_init handles init via state%soilwater
      ! hbot               = 0.0d0
      hbotab             = 0.0d0 
      hcomp              = 0.0d0 
      hdrain             = 0.0d0 
      ! [SS-SWC] hm1 retired — soilwater_init handles init via state%soilwater%hm1
      ! hm1                = 0.0d0
      ! [SS-CRP C-2.5] hroot retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! hroot              = 0.0d0
      hsublay            = 0.0d0
      ! [SS-BND B-2.7] hsurf retired — soilwater_init handles init via state%soilwater
      ! hsurf              = 0.0d0
      ! [SS-ATM A-2.6] igrai retired to state%atmosphere%intr%igrai
      ! [SS-SWC] ies0/iet0/iew0/iintc retired — soilwater_init handles init via state%soilwater
      ! ies0               = 0.0d0
      ! iet0               = 0.0d0
      ! iew0               = 0.0d0
      ! iintc              = 0.0d0
      infres             = 0.0d0 
      inpola             = 0.0d0 
      inpolb             = 0.0d0 
      ! [SS-SWC] inq retired — soilwater_init handles init via state%soilwater%inq
      ! inq                = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: inqdra/inqdra_in/inqdra_out removed (now in state%surfacewater).
      ! Initialized via surfacewater_init with allocate+assign=0.
      ! [SS-SWC] inqrot/inqssdi/ipondbeg/iprec/iqbot/iqtdo/iqtup/iqdo/iqup retired — soilwater_init handles init via state%soilwater
      ! inqrot             = 0.0d0
      ! inqssdi            = 0.0d0
      ! ipondbeg           = 0.0d0
      ! iprec              = 0.0d0
      ! iqbot              = 0.0d0
      ! iqtdo              = 0.0d0
      ! iqtup              = 0.0d0
      ! iqdo(1:numnod+1)   = 0.0d0
      ! iqup(1:numnod+1)   = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: iqdra removed (state%surfacewater%iqdra initialized in surfacewater_state_t).
      ! [SS-SWC] iqrot/iqssdi/iqredwet/iqreddry/iqredsol/iqredfrs/iruno/irunoCN/irunon retired — soilwater_init handles init via state%soilwater
      ! iqrot              = 0.0d0
      ! iqssdi             = 0.0d0
      ! iqredwet           = 0.0d0
      ! iqreddry           = 0.0d0
      ! iqredsol           = 0.0d0
      ! iqredfrs           = 0.0d0
      ! iruno              = 0.0d0
      ! irunoCN            = 0.0d0
      ! irunon             = 0.0d0
      issnowbeg          = 0.0d0 
      ! [SS-SWC] ithetabeg/k retired — soilwater_init handles init via state%soilwater
      ! ithetabeg          = 0.0d0
      ! k                  = 0.0d0
      khbot              = 0.0d0 
      khtop              = 0.0d0 
      ! [SS-SWC] kmean retired — soilwater_init handles init via state%soilwater%kmean
      ! kmean              = 0.0d0
      ! [GR-BH Task 36] ksatfit/ksatexm retired — seeded via state%soilwater in swap_mod.f90
      ! ksatfit            = 0.0d0
      ! ksatexm            = 0.0d0
      ksatthr            = 0.0d0 
      kvbot              = 0.0d0 
      kvtop              = 0.0d0 
      ! [GR-BH Task 37] L retired to state%drainage%L
      ! [SS-ATM A-2.6] ldwet retired to state%atmosphere%ldwet
      ! [SS-CRP C-2.5] mfluxtable retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! mfluxtable         = 0.0d0
      ! [SS-CRP C-2.5] mflux retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! mflux              = 0.0d0
      paramvg            = 0.0d0 
      ! [SS-SWC] pegwl/pond/pondini/pondm1 retired — soilwater_init handles init via state%soilwater
      ! pegwl              = 0.0d0
      ! pond               = 0.0d0
      ! pondini            = 0.0d0
      ! pondm1             = 0.0d0
      pondmx             = 0.0d0 
      pondmxtab          = 0.0d0 
      ! [SS-SWC] q retired — soilwater_init handles init via state%soilwater%q
      ! q                  = 0.0d0
      ! [SS-BND B-2.7] qbot, qbot_nonfrozen retired — soilwater_init handles init via state%soilwater
      ! qbot               = 0.0d0
      qbotab             = 0.0d0
      ! qbot_nonfrozen     = 0.0d0
      ! qdra               = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      ! qdrain             = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      qdraincomp         = 0.0d0 
      qdrtab             = 0.0d0 
      ! SS-SWST Phase 2 Task 11 C3: qdrtot removed (state%surfacewater%qdrtot initialized in surfacewater_state_t).
      qimmob             = 0.0d0 
      ! [SS-CRP C-2.5] qrosum/qred*sum/qrot retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! qrosum             = 0.0d0
      ! qredwetsum         = 0.0d0
      ! qreddrysum         = 0.0d0
      ! qredsolsum         = 0.0d0
      ! qredfrssum         = 0.0d0
      ! qrot               = 0.0d0
      qssdi              = 0.0d0
      ! dt_SSDI_event = 1.0 means "no SSDI sub-day event in progress";
      ! timecontrol.f90:423,620 read this ungated by flSSDI, so the
      ! default must reflect the no-event state even when SSDI is off.
      dt_SSDI_event      = 1.0d0
      ! [SS-BND B-2.7] qtop retired — soilwater_init handles init via state%soilwater
      ! qtop               = 0.0d0
      relsatthr          = 0.0d0 
      ! [SS-BND B-2.7] reva retired — soilwater_init handles init via state%soilwater
      ! reva               = 0.0d0
      ! [SS-HEAT] Task 9: rfcp retired to state%heat%rfcp
      ! rfcp               = 0.0d0
      rimlay             = 0.0d0 
      rsigni             = 0.0d0 
      rsoil              = 0.0d0 
      rsro               = 0.0d0 
      rsroexp            = 0.0d0
      Runoff_CN          = 0.0d0
      ! [SS-SWC] runon retired — soilwater_init handles init via state%soilwater%runon
      ! runon              = 0.0d0
      runonarr           = 0.0d0 
      ! [SS-BND B-2.7] runots retired — soilwater_init handles init via state%soilwater
      ! runots             = 0.0d0
      ! [SS-ATM A-2.6] saev retired to state%atmosphere%saev
      shape              = 0.0d0 
      sinamp             = 0.0d0 
      sinave             = 0.0d0 
      sinmax             = 0.0d0 
      ! [SS-ATM A-2.6] spev retired to state%atmosphere%spev
      sptab              = 0.0d0  
      sptablay           = 0.0d0 
      tau                = 0.0d0 
      ! [SS-SWC] theta/thetar/thetas/thetm1/thetsl/volact/volini/volm1 retired — soilwater_init handles init via state%soilwater
      ! theta              = 0.0d0
      ! thetar             = 0.0d0
      ! thetas             = 0.0d0
      ! thetm1             = 0.0d0
      ! thetsl             = 0.0d0
      ! volact             = 0.0d0
      ! volini             = 0.0d0
      ! volm1              = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: vtair removed (state%surfacewater%vtair initialized in surfacewater_state_t).
      ! [SS-SWC] wbalance retired — soilwater_init handles init via state%soilwater%wbalance
      ! wbalance           = 0.0d0
      ! wetper             = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      ! [GR-BH Task 35] z/ztopcp/zbotcp retired to state%mesh%z/ztopcp/zbotcp
      ! [SS-HEAT] Task 9: zfrostbot/zfrosttop retired to state%heat
      ! zfrostbot          = 0.0d0
      ! zfrosttop          = 0.0d0
      ! [GR-BH Task 37] zbotdr retired to state%drainage%zbotdr
      zi                 = 0.0d0
      zintf              = 0.0d0 
      ! ztopdislay         = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      ! [SS-TC] fldrain retired to state%timecontrol%flDrain (ADR 0041)
      ! [SS-SWC] fllowgwl retired — soilwater_init handles init via state%soilwater%fllowgwl
      ! fllowgwl           = .false.
      flrunon            = .false. 
      ! [SS-BND B-2.7] ftoph retired — soilwater_init handles init via state%soilwater
      ! ftoph              = .false.
!      drfil              = 16*' ' 
!      pathdrain          = 80*' '

! --- heat variables
      ! [SS-GR-FINAL D1] nheat zero-fill dropped — declaration retired 
      swbotbhea          = 0 
      swtopbhea          = 0 
      swcalt             = 0 
      swhea              = 0 
      ! [SS-GR-CROPRT A3] swtem zero-fill dropped — declaration retired
      tem                = 0 
      ddamp              = 0.0d0 
      ! [SS-HEAT] Task 9: fclay/forg/fquartz retired to state%heat
      ! fclay              = 0.0d0
      ! forg               = 0.0d0
      ! fquartz            = 0.0d0
      ! [GR-BH Task 36] orgmat/pclay/psand/psilt retired — seeded via state%soilwater in swap_mod.f90
      ! orgmat             = 0.0d0
      ! pclay              = 0.0d0
      ! psand              = 0.0d0
      ! psilt              = 0.0d0
      tampli             = 0.0d0 
      ! [SS-HEAT] Task 9: tebot retired to state%heat%tebot
      ! tebot              = 0.0d0
      tembtab            = 0.0d0
      temtoptab          = 0.0d0 
      tfroststa          = 0.0d0 
      tfrostend          = 0.0d0 
      timref             = 0.0d0 
      tmean              = 0.0d0 
      tsoil              = 0.0d0   ! [SS-HEAT] Task 9: config staging buffer — compute state is state%heat%tsoil
      ! [SS-HEAT] Task 9: tetop retired to state%heat%tetop
      ! tetop              = 0.0d0
      zh                 = 0.0d0 
      ! [SS-TC] fltemperature retired to state%timecontrol%flTemperature (ADR 0041)

! --- snow variables
      snw                = 0
      swsnow             = 0
      ! [SS-ATM A-2.6] cgsnow/cmelt/csnrai/csubl/gsnow/igsnow/isnrai/isubl/melt/
      !                snowinco/snrai/ssnow/slw/subl retired to state%atmosphere (cumu/intr/flat scalars)
      snowcoef           = 0.0d0
      TePrRain           = 0.0d0 
      TePrSnow           = 0.0d0 
      ! [SS-TC] flsnow retired to state%timecontrol%flSnow (ADR 0041)

! --- solute variables
      nconc              = 0
      ! [SS-GR-FINAL D1] sba zero-fill dropped — declaration retired 
      swbr               = 0 
      swsalinity         = 0
      swsolu             = 0 
      ! [SS-GR-CROPRT A1] AgeGwl1m zero-fill dropped — declaration retired (ADR 0032)
      bexp               = 0.0d0 
      cdrain             = 0.0d0 
      cirr               = 0.0d0 
      cml                = 0.0d0 
      cmsy               = 0.0d0 
!     cpond              = 0.0d0   ! moved to solute_state_t (ADR 0032)
      cpre               = 0.0d0
      cref               = 0.0d0
!     cseep              = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     csurf              = 0.0d0   ! moved to solute_state_t (ADR 0032)
      daquif             = 0.0d0 
      ddif               = 0.0d0 
      decpot             = 0.0d0 
      decsat             = 0.0d0 
!     dectot             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imdectot           = 0.0d0   ! moved to solute_state_t (ADR 0032)
      dtsolu             = 0.0d0 
      fdepth             = 0.0d0 
      frexp              = 0.0d0 
      gampar             = 0.0d0 
      ! [SS-GR-CROPRT A1] icAgeBot/icAgeDra/icAgeRot/icAgeSur zero-fills dropped — declarations retired (ADR 0032)
      isqbot             = 0.0d0 
      isqtop             = 0.0d0 
      kf                 = 0.0d0 
      kfsat              = 0.0d0 
      ldis               = 0.0d0 
      poros              = 0.0d0 
      rottot             = 0.0d0 
!     imrottot           = 0.0d0   ! moved to solute_state_t (ADR 0032)
      rtheta             = 0.0d0 
      salthead           = 0.0d0 
      saltmax            = 0.0d0 
      saltslope          = 0.0d0 
!     samcra             = 0.0d0   ! moved to solute_state_t (ADR 0032)
      samini             = 0.0d0
!     sampro             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     solbal             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqbot              = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imsqbot            = 0.0d0   ! moved to solute_state_t (ADR 0032)
      sqdra              = 0.0d0
!     imsqdra            = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqirrig            = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imsqirrig          = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqprec             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imsqprec           = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqrap              = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqsur              = 0.0d0   ! moved to solute_state_t (ADR 0032)
      tscf               = 0.0d0 
      zc                 = 0.0d0 
      ! [SS-TC] flsolute retired to state%timecontrol%flSolute (ADR 0041)
      ! [SS-GR-CROPRT A1] flAgeTracer zero-fill dropped — declaration retired (ADR 0032)

! --- macropore retirement [MACRO-RETIRE 2026-05-12] — ADR 0040
!     Only the retired-zero placeholders kept in variables.f90 are
!     initialized here. See legacy/swap-4.2.0 for the original
!     macropore initialization that populated ~80 macropore globals.
      ! [SS-GR-FINAL D1] SwSoilShr zero-fill dropped — declaration retired
      ! [SS-GR-FINAL D1] ThetCrMp zero-fill dropped — declaration retired
      Z_Tp                 = 0.0d0
      CritUndSatVol        = 0.0d0
      ! [GR-BH Task 36] ArMpSs retired — made local in soilhydraulics/solute/agetracer (ADR 0040 complete)
      ! ArMpSs               = 0.0d0
      ArMpTp               = 0.0d0
      cQMpLatSs            = 0.0d0
      cQMpOutDrRap         = 0.0d0
      ! [SS-GR-FINAL D1] dFdhMp zero-fill dropped — declaration retired
      iQMpOutDrRap         = 0.0d0
      ! [SS-GR-CROPRT A2] iQInTopLatDm1/2/iQInTopVrtDm1/2 zero-fills dropped — declarations retired (ADR 0040)
      ! [SS-GR-CROPRT A2] IWaSrDm1/2Beg/WaSrDm1/2/WaSrDm1/2Ini zero-fills dropped — declarations retired (ADR 0040)
      DiPoCp               = 0.0d0
      IAvFrMpWlWtDm1       = 0.0d0
      IAvFrMpWlWtDm2       = 0.0d0
      iQExcMtxDm1Cp        = 0.0d0
      iQExcMtxDm2Cp        = 0.0d0
      iQOutDrRapCp         = 0.0d0
      VlMpStDm1            = 0.0d0
      VlMpStDm2            = 0.0d0
      ! [SS-GR-CROPRT A2] IcTopMP zero-fill dropped — declaration retired (ADR 0040)
      ! [SS-GR-FINAL D1] IDecMpRat zero-fill dropped — declaration retired
      QExcMpMtx            = 0.0d0
      QMaPo                = 0.0d0
      QRapDra              = 0.0d0
      ! Drainage globals (NOT macropore — kept here next to retired set)
      NumLevRapDra         = 0
      ! [SS-GR-FINAL D1] RapDraReaExp/RapDraResRef zero-fills dropped — declarations retired
      FlDecMpRat           = .false.
      ! [SS-GR-CROPRT A2] flmacropore zero-fill dropped — declaration retired (ADR 0040)

! --- surface water variables
      ! [SS-GR-FINAL D1] swswb/swdrf zero-fills dropped — declarations retired
      swsrf                = 0
      swallo               = 0
      swdtyp               = 0
      ! [GR-BH Task 37] swnrsrf retired to state%drainage%swnrsrf
      swqhr                = 0
      swsec                = 0
      nrpri                = 0
      ! [SS-GR-FINAL D1] nrsec zero-fill dropped — declaration retired
      nmper                = 0
      swman                = 0
      ! [GR-BH Task 37] SwTopnrsrf retired to state%drainage%swtopnrsrf
      ! [SS-GR-FINAL D1] nqh/drf/swb zero-fills dropped — declarations retired
      nphase               = 0
      nodhd                = 0
      ! SS-SWST Phase 2 Task 11 C3: numadj/imper removed (state%surfacewater owns them).
      intwl                = 0
      ! numadj=0, imper=1 set by surfacewater_state_t defaults.
      nowltab              = 0
      widthr                = 0.0d0
      taludr                = 0.0d0
      rdrain                = 0.0d0
      rsurfdeep             = 0.0d0
      rsurfshallow          = 0.0d0
      rinfi                 = 0.0d0
      rentry                = 0.0d0
      rexit                 = 0.0d0
      gwlinf                = 0.0d0
      wlptab                = 0.0d0
      impend                = 0.0d0
      wldip                 = 0.0d0
      wscap                 = 0.0d0
      hbweir                = 0.0d0
      osswlm                = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: wlstar removed (state%surfacewater%wlstar, default=0).
      wlp                   = 0.0d0
      alphaw                = 0.0d0
      betaw                 = 0.0d0
      dropr                 = 0.0d0
      hdepth                = 0.0d0
      gwlcrit               = 0.0d0
      hcrit                 = 0.0d0
      vcrit                 = 0.0d0
      hqhtab                = 0.0d0
      qqhtab                = 0.0d0
      wlsman                = 0.0d0
      wlstab                = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: sttab/swstini/swst/wlsbak removed (state%surfacewater owns them).
      ! sttab=0, swstini=0, swst=0, wlsbak=0 are set by surfacewater_state_t defaults.
      cofintfl              = 0.0d0
      expintfl              = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: cqdrd/cwsupp/cwout/wls removed (state%surfacewater owns them).
      ! All default to 0 in surfacewater_state_t.
      owltab                = 0.0d0
      ! SS-SWST Phase 2 Task 11 C3: hwlman/wlsold removed (state%surfacewater owns them, default=0).
      ! qdrd                  = 0.0d0  ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      ! [SS-TC] flsurfacewater retired to state%timecontrol%flSurfaceWater (ADR 0041)
      ! SS-SWST Phase 2 Task 11 C3: overfl removed (state%surfacewater%overfl, default=.false.).

      return
      end subroutine
    
      Subroutine InitializeCrop
! --- Initialize all crop specific variables in subroutine Variables
      use variables
      implicit none 
      
! --- Initialize all crop variables 

!     date calendar starts is date of workability or crop start
      flCropPrep         = .false.
      flCropSow          = .false.
      flCropGerm         = .false.
      
      PrepDelay          =  0
      SowDelay           =  0

      tsumgerm           =  0.d0

      flCropHarvest      = .false.
      flHarvestDay       = .false.
      ! [SS-CRP C-2.5] flWrtNonox retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! flWrtNonox         = .false.

      daycrop            = 0 
      daygrowth          = 0
      daygrowthpot       = 0
      idev               = 0 
      idsl               = 0 
      noddrz             = 0 
      swcf               = 0 
      swgc               = 0 
      swoxygen           = 0 
      adcrh              = 0.0d0 
      adcrl              = 0.0d0 
      albedo             = 0.0d0 
      alphacrit          = 0.0d0 
      amaxtb             = 0.0d0 
      atmin7             = 0.0d0
      avevaptb           = 0.0d0 
      avprectb           = 0.0d0 
      cf                 = 0.0d0 
      cfeictb            = 0.0d0 
      cftb               = 0.0d0 
      ch                 = 0.0d0 
      chtb               = 0.0d0 
      gctb               = 0.0d0 
      cofab              = 0.0d0 
      cumdens            = 0.0d0 
      cuptgraz           = 0.0d0 
      cuptgrazpot        = 0.0d0 
      cvl                = 0.0d0 
      cvo                = 0.0d0 
      cvr                = 0.0d0 
      cvs                = 0.0d0 
      cwdm               = 0.0d0 
      cwdmpot            = 0.0d0 
      dlc                = 0.0d0 
      dlo                = 0.0d0 
      dtsmtb             = 0.0d0
      dvsend             = 0.0d0
      eff                = 0.0d0 
      fbl                = 0.0d0 
      fbltb              = 0.0d0 
      fltb               = 0.0d0 
      fotb               = 0.0d0 
      frtb               = 0.0d0 
      fstb               = 0.0d0 
      gasstpot           = 0.0d0 
      gc                 = 0.0d0 
      hlim1              = 0.0d0 
      hlim2l             = 0.0d0 
      hlim2u             = 0.0d0 
      hlim3h             = 0.0d0 
      hlim3l             = 0.0d0 
      hlim4              = 0.0d0 
      kdif               = 0.0d0 
      kdir               = 0.0d0
      laiem              = 0.0d0
      lv                 = 0.0d0
      lvage              = 0.0d0
      lvagepot           = 0.0d0
      lvpot              = 0.0d0
      mrestpot           = 0.0d0 
      mrftb              = 0.0d0
      perdl              = 0.0d0 
      pfreetb            = 0.0d0 
      pstemtb            = 0.0d0 
      q10                = 0.0d0 
      rdrrtb             = 0.0d0 
      rdrstb             = 0.0d0 
      reltr              = 0.0d0 
      rfsetb             = 0.0d0 
      rgrlai             = 0.0d0 
      rml                = 0.0d0 
      rmo                = 0.0d0 
      rmr                = 0.0d0 
      rms                = 0.0d0 
      rootcoefa          = 0.0d0 
      rootradius         = 0.0d0 
      ! [SS-CRP C-2.5] rootrho retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! rootrho            = 0.0d0
      rsc                = 0.0d0
      rsw                = 0.0d0 

      scanopytb          = 0.0d0 
      siccaptb           = 0.0d0 
      sla                = 0.0d0
      slapot             = 0.0d0
      slatb              = 0.0d0 
      spa                = 0.0d0 
      span               = 0.0d0 
      ssa                = 0.0d0
      tagppot            = 0.0d0
      tagpt              = 0.0d0 
      tagptpot           = 0.0d0 
      tbase              = 0.0d0 
      tdwi               = 0.0d0 
      tmnftb             = 0.0d0 
      tmpftb             = 0.0d0 
      tsumam             = 0.0d0 
      tsumea             = 0.0d0 
      vernrtb            = 0.0d0 

! --- only for bulb crops (tulips etc..)
      wbl                = 0.0d0 
      wblpot             = 0.0d0
      dwbl               = 0.0d0
      dwblpot            = 0.0d0
      plwt               = 0.0d0
      
      wiltpoint          = 0.0d0 
      wrtb               = 0.0d0
      zsow               = 0.0d0
      wrtpot             = 0.0d0
      dwst               = 0.0d0
      dwrt               = 0.0d0
      dwso               = 0.0d0
      
      dwlvCrop           = 0.0d0
      dwlvSoil           = 0.0d0
      HarLosOrm_tot      = 0.0d0
      
! --- interception
      ! [SS-ATM A-2.6] sicact retired to state%atmosphere%sicact
      
! --- rooting
      rdctb              = 0.0d0 
      rdtb               = 0.0d0
      rlwtb              = 0.0d0
      rdi                = 0.0d0 
      rri                = 0.0d0       
      rdc                = 0.0d0 
      wrtmax             = 0.0d0

! --- harvest grassland
      dateharvest        = 0.0d0
      dmmowtb            = 0.0d0
      dmgrztb            = 0.0d0
      lsda               = 0.0d0
      DelayRegrowthTab   = 0.0d0
      DaysGrazingTab     = 0.0d0
      UptGrazingTab      = 0.0d0
      LossGrazingTab     = 0.0d0
      lossmowtab         = 0.0d0
      lossgrztab         = 0.0d0
      
      flanthesis         = .false.
      flgrazing          = .false.
      flgrazingpot       = .false.
      flharvest          = .false.
      flharvestpot       = .false. 

      return
      end subroutine
