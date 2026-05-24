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
      ! [GR-ATM 2026-05-23] logf retired — swap_log owns the log-file unit
      ! [SS-GR-FINAL D1] ex_tlast zero-fill dropped — declaration retired
      outdat             = 0.0d0
      outdatint          = 0.0d0
      flSwapShared       = .false.

! --- crop variables (not crop specific, such as crop calendar)
      crp                = 0 
      icrop              = 0 
      ! [GR-ATM 2026-05-23] croptype retired — see state%crop%common%croptype
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
      swetsine           = 0
      ! [SS-ATM A-2.6] aintcdt/atmdem/caintc/cevap/cgrai/cnrai/cpeva/cptra retired to state%atmosphere
      ! [GR-ATM 2026-05-23] daylp retired — see state%atmosphere%daylp

! --- only for bulb crops (tulips etc..)
      ! [SS-ATM A-2.6] empreva/fprecnosnow/grai/graidt/ievap/inrai/ipeva/iptra retired to state%atmosphere
      epot               = 0.0d0
      ! [GR-ATM 2026-05-23] cfevappond retired — see state%atmosphere%cfevappond
      ! [GR-ATM 2026-05-23] finterception retired — see state%atmosphere%finterception
      grain              = 0.0d0
      ! [SS-TC] metperiod retired to state%timecontrol%metperiod (ADR 0041)
      ! [SS-ATM A-2.6] nraida/nraidt/peva/pevaday/ptra/ptraday retired to state%atmosphere
      nrain              = 0.0d0
      ! [GR-ATM 2026-05-23] rad retired — see state%atmosphere%rad
      ! [GR-CROP C12] rainamount/rainfluxarray/raintimearray retired from variables — state%atmosphere%X
      ! [GR-ATM C8] rh retired — state%atmosphere%rh written by meteoday
      tadw               = 0.0d0
      tadwpot            = 0.0d0
      ! [GR-ATM 2026-05-23] tav retired — see state%atmosphere%Tav
      ! [GR-ATM C8] tavd retired — state%atmosphere%tavd written by meteoday
      ! [GR-ATM 2026-05-23] tmn retired — see state%atmosphere%tmn
      tmnr               = 0.0d0 
      ! [GR-ATM 2026-05-23] tmx retired — see state%atmosphere%tmx 
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
      swirfix            = 0 
      cirrs              = 0.0d0 
      dcrit              = 0.0d0 
      ditab              = 0.0d0 
      dwatab             = 0.0d0 
      fidtab             = 0.0d0 
      hcritab            = 0.0d0 
      ! [SS-SWC] igird/inird retired — soilwater_init handles init via state%soilwater
      ! igird              = 0.0d0
      ! inird              = 0.0d0
      irconc             = 0.0d0 
      irdate             = 0.0d0 
      irdepth            = 0.0d0 
      ! [GR-ATM 2026-05-23] nird retired — see state%atmosphere%nird
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
      ! [SS-BMI2 Task 5] MaxIt retired 2026-05-13 — moved to state%timecontrol%MaxIt
      Itnumb             = 0
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
      ! [GR-SOIL 2026-05-24] botcom retired — see state%mesh%botcom
      dra                = 0 
      ! [GR-DRA 2026-05-23] dramet retired — see state%drainage%dramet
      inc                = 0 
      ! [SS-SWC] indeks retired — soilwater_init handles init via state%soilwater%indeks
      ! indeks             = 0
      ! [GR-DRA 2026-05-23] ipos retired — see state%drainage%ipos
      ! [GR-SOIL 2026-05-24] isoillay retired — read inline from config%soil%isoillay
      ! [GR-BH Task 35] layer retired to state%mesh%layer
      ! [SS-BMI2 Task 5] msteps retired 2026-05-13 — moved to state%timecontrol%msteps
      ! [GR-SOIL 2026-05-24] ncomp retired — read inline from config%soil%ncomp
      nhead              = 0 
      ! [SS-SWC] nodgwl retired — soilwater_init handles init via state%soilwater%nodgwl
      ! nodgwl             = 0
      ! [GR-SOIL 2026-05-24] nod1lay retired — see state%mesh%nod1lay
      ! [SS-HEAT] Task 9: nodfrostbot retired to state%heat%nodfrostbot
      ! nodfrostbot        = 0
      ! [SS-SWC] npegwl retired — soilwater_init handles init via state%soilwater%npegwl
      ! npegwl             = 0
      ! [GR-BH Task 37] nrlevs retired to state%drainage%nrlevs
      ! [GR-SOIL 2026-05-24] nsublay retired — derived inline from config%soil%sublay
      numbit             = 0 
      ! [GR-SOIL 2026-05-24] numlay retired — see state%mesh%numlay
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
      ! [GR-SOL 2026-05-24] swbotbc retired — see state%solute%swbotbc
      swbotb3Impl        = 0 
      ! [GR-BND 2026-05-23] SwBotb3ResVert retired — see state%soilwater%swbotb3resvert
      swcfbs             = 0 
      swdiscrvert        = 0 
      ! [GR-BH Task 37] swdivd retired to state%drainage%swdivd
      ! [GR-BH Task 37] swdivdinf retired to state%drainage%swdivdinf
      ! [GR-DRA 2026-05-23] swdislay/swtopdislay retired — see state%drainage%X
      swfrost            = 0
      swhyst             = 0 
      ! [GR-SOL 2026-05-24] swinco retired — see state%soilwater%swinco
      ! [SS-GR-FINAL D5] swoutputmodflow zero-fill dropped — declaration retired
      ! [GR-BND 2026-05-23] swpondmx retired — see state%surfacewater%swpondmx
      swqhbot            = 0 
      ! [GR-ATM 2026-05-23] swredu retired — see state%atmosphere%swredu
      ! [SS-GR-CROPRT C1] swend zero-fill dropped — global + state field retired (ADR 0009: always 0)
      ! [SS-GR-FINAL D1] swvap/vap/wba zero-fills dropped — declarations retired
      aqamp              = 0.0d0 
      aqave              = 0.0d0 
      aqper              = 0.0d0 
      aqtmax             = 0.0d0 
      ! [GR-DRA 2026-05-23] basegw retired — see state%drainage%basegw
      ! [GR-SOL 2026-05-24] bdens retired — see state%soilwater%bdens
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
      ! [GR-ATM 2026-05-23] cofred retired — see state%atmosphere%cofred
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
      ! [GR-SOL 2026-05-24] cseeptab retired — see state%solute%cseeptab
      ! [SS-BND B-2.7] deepgw retired — soilwater_init handles init via state%soilwater
      ! deepgw             = 0.0d0
      ! [SS-SWC] dimoca retired — soilwater_init handles init via state%soilwater%dimoca
      ! dimoca             = 0.0d0
      ! [GR-BH Task 35] disnod retired to state%mesh%disnod
      ! drainl             = 0.0d0  ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      drares             = 0.0d0 
      ! [GR-BH Task 35] dz retired to state%mesh%dz
      dznew              = 0.0d0 
      ! [GR-DRA 2026-05-23] entres retired — see state%drainage%entres
      ! [SS-SWC] evp retired — soilwater_init handles init via state%soilwater%evp
      ! evp                = 0.0d0
      ! [GR-BH Task 37] FacDpthInf retired to state%drainage%FacDpthInf
      ! [GR-DRA 2026-05-23] ftopdislay retired — see state%drainage%ftopdislay
      ! [GR-DRA 2026-05-23] geofac retired — see state%drainage%geofac
      ! [SS-SWC] gwl retired — soilwater_init handles init via state%soilwater%gwl
      ! gwl                = 0.0d0
      gwli               = 0.0d0 
      ! [SS-BND B-2.7] gwlinp retired — soilwater_init handles init via state%soilwater
      ! gwlinp             = 0.0d0
      ! [SS-SWC] gwlm1 retired — soilwater_init handles init via state%soilwater%gwlm1
      ! gwlm1              = 0.0d0
      ! [GR-BND 2026-05-23] gwltab retired — soilwater_init zeros state%soilwater%gwltab
      ! [SS-SWC] h retired — soilwater_init handles init via state%soilwater%h
      ! h                  = 0.0d0
      h_enpr             = 0.0d0 
      ! [GR-BND 2026-05-23] haqtab retired — soilwater_init zeros state%soilwater%haqtab
      ! [SS-SWC] hatm retired — soilwater_init handles init via state%soilwater%hatm
      ! hatm               = 0.0d0
      ! [SS-BND B-2.7] hbot retired — soilwater_init handles init via state%soilwater
      ! hbot               = 0.0d0
      ! [GR-BND 2026-05-23] hbotab retired — soilwater_init zeros state%soilwater%hbotab
      ! [GR-SOIL 2026-05-24] hcomp retired — derived inline from config%soil%hsublay/ncomp
      hdrain             = 0.0d0 
      ! [SS-SWC] hm1 retired — soilwater_init handles init via state%soilwater%hm1
      ! hm1                = 0.0d0
      ! [SS-CRP C-2.5] hroot retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! hroot              = 0.0d0
      ! [GR-SOIL 2026-05-24] hsublay retired — read inline from config%soil%hsublay
      ! [SS-BND B-2.7] hsurf retired — soilwater_init handles init via state%soilwater
      ! hsurf              = 0.0d0
      ! [SS-ATM A-2.6] igrai retired to state%atmosphere%intr%igrai
      ! [SS-SWC] ies0/iet0/iew0/iintc retired — soilwater_init handles init via state%soilwater
      ! ies0               = 0.0d0
      ! iet0               = 0.0d0
      ! iew0               = 0.0d0
      ! iintc              = 0.0d0
      infres             = 0.0d0 
      ! [GR-SOL 2026-05-24] inpola/inpolb retired — see state%mesh%{inpola,inpolb}
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
      ! [GR-ATM 2026-05-23] issnowbeg retired — see state%atmosphere%ISsnowBeg
      ! [SS-SWC] ithetabeg/k retired — soilwater_init handles init via state%soilwater
      ! ithetabeg          = 0.0d0
      ! k                  = 0.0d0
      ! [GR-DRA 2026-05-23] khbot retired — see state%drainage%khbot
      ! [GR-DRA 2026-05-23] khtop retired — see state%drainage%khtop
      ! [SS-SWC] kmean retired — soilwater_init handles init via state%soilwater%kmean
      ! kmean              = 0.0d0
      ! [GR-BH Task 36] ksatfit/ksatexm retired — seeded via state%soilwater in swap_mod.f90
      ! ksatfit            = 0.0d0
      ! ksatexm            = 0.0d0
      ksatthr            = 0.0d0 
      ! [GR-DRA 2026-05-23] kvbot retired — see state%drainage%kvbot
      ! [GR-DRA 2026-05-23] kvtop retired — see state%drainage%kvtop
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
      ! [GR-BND 2026-05-23] pondmxtab retired — surfacewater_state_init zeros state%surfacewater%pondmxtab
      ! [SS-SWC] q retired — soilwater_init handles init via state%soilwater%q
      ! q                  = 0.0d0
      ! [SS-BND B-2.7] qbot, qbot_nonfrozen retired — soilwater_init handles init via state%soilwater
      ! qbot               = 0.0d0
      ! [GR-BND 2026-05-23] qbotab retired — soilwater_init zeros state%soilwater%qbotab
      ! qbot_nonfrozen     = 0.0d0
      ! qdra               = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      ! qdrain             = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      qdraincomp         = 0.0d0 
      ! [GR-DRA 2026-05-23] qdrtab retired — see state%drainage%qdrtab
      ! SS-SWST Phase 2 Task 11 C3: qdrtot removed (state%surfacewater%qdrtot initialized in surfacewater_state_t).
      ! [GR-SOIL 2026-05-24] qimmob retired — fingered-flow flux retired-zero inline in waterbalance.f90.
      ! [SS-CRP C-2.5] qrosum/qred*sum/qrot retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! qrosum             = 0.0d0
      ! qredwetsum         = 0.0d0
      ! qreddrysum         = 0.0d0
      ! qredsolsum         = 0.0d0
      ! qredfrssum         = 0.0d0
      ! qrot               = 0.0d0
      ! [GR-SOIL 2026-05-24] qssdi init handled by soilwater_init (state field).
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
      ! [GR-ATM 2026-05-23] rsigni retired — see state%atmosphere%rsigni
      ! [GR-ATM 2026-05-23] rsoil retired — see state%atmosphere%rsoil
      ! [GR-ATM 2026-05-23] Runoff_CN retired — see state%atmosphere%Runoff_CN
      ! [SS-SWC] runon retired — soilwater_init handles init via state%soilwater%runon
      ! runon              = 0.0d0
      ! [GR-BND 2026-05-23] runonarr retired — soilwater_init zeros state%soilwater%runonarr
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
      ! [GR-DRA 2026-05-23] zintf retired — see state%drainage%zintf 
      ! ztopdislay         = 0.0d0   ! Moved to drainage_state_t — drainage_init handles (ADR 0031)
      ! [SS-TC] fldrain retired to state%timecontrol%flDrain (ADR 0041)
      ! [SS-SWC] fllowgwl retired — soilwater_init handles init via state%soilwater%fllowgwl
      ! fllowgwl           = .false.
      ! [GR-BND 2026-05-23] flrunon retired — see state%soilwater%flrunon
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
      ! [SS-TC] flsnow retired to state%timecontrol%flSnow (ADR 0041)

! --- solute variables
      ! [GR-SOL 2026-05-24] nconc retired — see state%solute%nconc
      ! [SS-GR-FINAL D1] sba zero-fill dropped — declaration retired 
      ! [GR-SOL 2026-05-24] swbr retired — see state%solute%swbr
      swsolu             = 0 
      ! [SS-GR-CROPRT A1] AgeGwl1m zero-fill dropped — declaration retired (ADR 0032)
      ! [GR-SOL 2026-05-24] bexp retired — see state%solute%bexp
      ! [GR-SOL 2026-05-24] cdrain retired — see state%solute%cdrain
      ! [GR-SOL 2026-05-24] cirr retired — see state%solute%cirr
      ! [GR-SOL 2026-05-24] cml retired — see state%solute%cml_init
      cmsy               = 0.0d0 
!     cpond              = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] cpre retired — see state%solute%cpre
      ! [GR-SOL 2026-05-24] cref retired — see state%solute%cref
!     cseep              = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     csurf              = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] daquif retired — see state%solute%daquif
      ! [GR-SOL 2026-05-24] ddif retired — see state%solute%ddif
      ! [GR-SOL 2026-05-24] decpot retired — see state%solute%decpot
      ! [GR-SOL 2026-05-24] decsat retired — see state%solute%decsat
!     dectot             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imdectot           = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] dtsolu retired — see state%solute%dtsolu
      ! [GR-SOL 2026-05-24] fdepth retired — see state%solute%fdepth
      ! [GR-SOL 2026-05-24] frexp retired — see state%solute%frexp
      ! [GR-SOL 2026-05-24] gampar retired — see state%solute%gampar
      ! [SS-GR-CROPRT A1] icAgeBot/icAgeDra/icAgeRot/icAgeSur zero-fills dropped — declarations retired (ADR 0032)
      ! [GR-SOL 2026-05-24] isqbot retired — see state%solute%isqbot
      ! [GR-SOL 2026-05-24] isqtop retired — see state%solute%isqtop
      ! [GR-SOL 2026-05-24] kf retired — see state%solute%kf
      ! [GR-SOL 2026-05-24] kfsat retired — see state%solute%kfsat
      ! [GR-SOL 2026-05-24] ldis retired — see state%solute%ldis
      ! [GR-SOL 2026-05-24] poros retired — see state%solute%poros
      ! [GR-SOL 2026-05-24] rottot retired — see state%solute%rottot
!     imrottot           = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] rtheta retired — see state%solute%rtheta
!     samcra             = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] samini retired — see state%solute%samini
!     sampro             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     solbal             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqbot              = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imsqbot            = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] sqdra retired — see state%solute%sqdra
!     imsqdra            = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqirrig            = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imsqirrig          = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqprec             = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     imsqprec           = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqrap              = 0.0d0   ! moved to solute_state_t (ADR 0032)
!     sqsur              = 0.0d0   ! moved to solute_state_t (ADR 0032)
      ! [GR-SOL 2026-05-24] tscf retired — see state%solute%tscf
      ! [GR-SOL 2026-05-24] zc retired — see state%solute%zc_init
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
      ! [GR-SOIL 2026-05-24] QExcMpMtx + QMaPo retired — ADR 0040 macropore terms retired-zero inline in waterbalance.f90.
      ! [GR-DRA 2026-05-23] QRapDra retired — see state%drainage%QRapDra
      ! Drainage globals (NOT macropore — kept here next to retired set)
      ! [GR-DRA 2026-05-23] NumLevRapDra retired — see state%drainage%NumLevRapDra
      ! [SS-GR-FINAL D1] RapDraReaExp/RapDraResRef zero-fills dropped — declarations retired
      FlDecMpRat           = .false.
      ! [SS-GR-CROPRT A2] flmacropore zero-fill dropped — declaration retired (ADR 0040)

! --- surface water variables
      ! [SS-GR-FINAL D1] swswb/swdrf zero-fills dropped — declarations retired
      ! [GR-DRA 2026-05-23] swallo/swdtyp retired — see state%drainage%{swallo,swdtyp}
      ! [GR-BH Task 37] swnrsrf retired to state%drainage%swnrsrf
      ! [GR-DRA 2026-05-23] swqhr retired — see state%surfacewater%swqhr
      ! [GR-DRA 2026-05-23] nrpri retired — see state%surfacewater%nrpri
      ! [SS-GR-FINAL D1] nrsec zero-fill dropped — declaration retired
      ! [GR-DRA 2026-05-23] nmper retired — see state%surfacewater%nmper
      ! [GR-DRA 2026-05-23] swman retired — see state%surfacewater%swman
      ! [GR-BH Task 37] SwTopnrsrf retired to state%drainage%swtopnrsrf
      ! [SS-GR-FINAL D1] nqh/drf/swb zero-fills dropped — declarations retired
      ! [GR-DRA 2026-05-23] nphase/nodhd retired — see state%surfacewater%{nphase,nodhd}
      ! SS-SWST Phase 2 Task 11 C3: numadj/imper removed (state%surfacewater owns them).
      ! [GR-DRA 2026-05-23] intwl retired — see state%surfacewater%intwl
      ! numadj=0, imper=1 set by surfacewater_state_t defaults.
      ! [GR-DRA 2026-05-23] nowltab retired — see state%drainage%nowltab
      widthr                = 0.0d0
      taludr                = 0.0d0
      rdrain                = 0.0d0
      rsurfdeep             = 0.0d0
      rsurfshallow          = 0.0d0
      rinfi                 = 0.0d0
      rentry                = 0.0d0
      rexit                 = 0.0d0
      gwlinf                = 0.0d0
      ! [GR-DRA 2026-05-23] wlptab retired — see state%surfacewater%wlptab
      ! [GR-DRA 2026-05-23] impend/wldip/wscap retired — see state%surfacewater%X
      ! [GR-DRA 2026-05-23] hbweir retired — see state%surfacewater%hbweir
      ! [GR-DRA 2026-05-23] osswlm retired — see state%surfacewater%osswlm
      ! SS-SWST Phase 2 Task 11 C3: wlstar removed (state%surfacewater%wlstar, default=0).
      wlp                   = 0.0d0
      ! [GR-DRA 2026-05-23] alphaw/betaw retired — see state%surfacewater%{alphaw,betaw}
      ! [GR-DRA 2026-05-23] dropr retired — see state%surfacewater%dropr
      hdepth                = 0.0d0
      ! [GR-DRA 2026-05-23] gwlcrit/hcrit/vcrit retired — see state%surfacewater%X
      hqhtab                = 0.0d0
      ! [GR-DRA 2026-05-23] qqhtab retired — see state%surfacewater%qqhtab
      ! [GR-DRA 2026-05-23] wlsman retired — see state%surfacewater%wlsman
      ! [GR-DRA 2026-05-23] wlstab retired — see state%surfacewater%wlstab
      ! SS-SWST Phase 2 Task 11 C3: sttab/swstini/swst/wlsbak removed (state%surfacewater owns them).
      ! sttab=0, swstini=0, swst=0, wlsbak=0 are set by surfacewater_state_t defaults.
      ! [GR-DRA 2026-05-23] cofintfl/expintfl retired — see state%drainage%{cofintfl,expintfl}
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
      idsl               = 0 
      noddrz             = 0 
      atmin7             = 0.0d0
      avevaptb           = 0.0d0 
      avprectb           = 0.0d0 
      dlc                = 0.0d0 
      dlo                = 0.0d0 
      fbl                = 0.0d0 
      gasstpot           = 0.0d0 
      ! [GR-ATM 2026-05-23] gc retired — see state%crop%common%gc
      mrestpot           = 0.0d0 
      mrftb              = 0.0d0
      pfreetb            = 0.0d0 
      pstemtb            = 0.0d0 
      reltr              = 0.0d0 
      rootcoefa          = 0.0d0 
      rootradius         = 0.0d0 
      ! [SS-CRP C-2.5] rootrho retired — soilwater_init handles init via state%soilwater (ADR 0036)
      ! rootrho            = 0.0d0

      scanopytb          = 0.0d0 
      ! [GR-ATM 2026-05-23] siccaptb retired — see state%crop%common%siccaptb
      vernrtb            = 0.0d0 

! --- only for bulb crops (tulips etc..)
      wiltpoint          = 0.0d0 
      wrtb               = 0.0d0
      zsow               = 0.0d0
! --- interception
      ! [SS-ATM A-2.6] sicact retired to state%atmosphere%sicact
      
! --- rooting

! --- harvest grassland
      
      flanthesis         = .false.

      return
      end subroutine
