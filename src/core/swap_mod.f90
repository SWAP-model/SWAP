!> @file swap_mod.f90
!! SS-DRV Phase 1: module form of the legacy `subroutine swap`.
!! Three named lifecycle procedures replace the (iCaller, iTask) dispatch.
!! Time loop lives in the caller. State and config are threaded explicitly.
!!
!! Note: heavy compute `use` statements are scoped to each procedure rather
!! than the module header. The unit-test build only provides state/config
!! modules, and hoisting all uses would require the full compute module chain
!! at compile time even for the lightweight smoke tests.
module swap_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private
   public :: swap_init, swap_run_step, swap_close, swap_init_from_loaded_config

contains

   subroutine swap_init(config_file, state, config)
      use load_swap_config_mod, only: load_swap_config
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config  ! target: crop_config_global => config%crop (set inside body, Task 3)

!  Phase 4f strangler-fig: read time-independent input via the TOML
!  pipeline + config_to_variables adapter. The legacy readswap() entry
!  point is no longer called from the runtime path (see ADR 0007); it
!  survives in src/io/readswap.f90 only as a parity-test fixture. The
!  binary expects swap.toml in the current directory; abort_if_fatal
!  terminates with a clear summary if the file is absent or fails
!  validate/finalize.
!
!  Remaining strangler-fig debt: ~10 individual HACK Phase 4f-extend
!  slots in config_to_variables.f90 for legacy globals not yet covered
!  by typed schema slots (SWREDU, RSIGNI, CFEVAPPOND, iHWCKmodel, RDS,
!  ksatexm path, etc.). Each slot is a small typed-config extension +
!  adapter wiring. The deeper follow-on (per ADR 0016) is the config-
!  passing refactor — eliminate variables-module mutation by passing
!  typed config + state to compute subs explicitly.
   block
      use error_mod, only: error_collection_t
      type(error_collection_t) :: errors
      call load_swap_config(config_file, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
      call errors%abort_if_fatal()
   end block

   call swap_init_from_loaded_config(state, config)

   end subroutine swap_init

   subroutine swap_init_from_loaded_config(state, config)
      use variables, only : flswapshared, flcropnut, swfrost, &
                            ! [SS-GR-CROPRT A1] flagetracer dropped — retired (ADR 0032; always .false.)
                            swusecn, flcropcalendar, &
                            flharvestday, flcropoutput, swcrp, project, &
                            ! [SS-GR-CROPRT C1] swend dropped — global + state field retired (ADR 0009; always 0)
                            flCropHarvest, &   ! [SS-GR-CROPRT A5] initial zero mirror
                            ! [GR-FINAL C3] flirrigationoutput dropped: W-global (0 consumers; ADR 0009 deleted IrrigationOutput)
                            flTillage, flSSDI, &
                            numlay, &
                            owltab, nowltab, &
                            ! [GR-CROP-DVS] meteo arrays/scalars + CN runoff tables retired
                            !   (arad,atmn,atmx,ahum,awin,arai,aetr,wet,atav,epot,tpot,grain,nrain,
                            !    tav,daynrfirst,daynrlast,atmin7,nofd,isua,avevaptb,avprectb,
                            !    pfreetb,pstemtb,scanopytb,CNdry,CNwet,ThetaRef,Runoff_CN,
                            !    wc_cor,wc10,iCNtab,CNtimTAB,CNrefTAB) → state%atmosphere

                            ! [GR-ATM C8] out_tmn/tmx/hum/win/etr/wet/rad retired from import (state written by meteoday)
                            swcfbs, flCropEmergence, &  ! lai/swcf retired
                            ! [SS-GR-CROP A14] crop_common legacy globals; dvs/tsum retired
                            daycrop, icrop, &
                            ! rd/rdpot/rdm/rdi/rri/rdc/HarLosOrm_tot/cuptgraz/cuptgrazpot retired
                            ! [SS-GR-CROP A15] crop_wofost/grass/fixed legacy globals; biomass+dw*+plossdm/lossdm retired
                            swbulb, &
                            ! mowrest/seqgrazmow/seqgrazmowpot/dateharvest retired
                            swpotrelmf  ! cropstart/end act/pot retired; cftb/chtb/cfeic/cfeictb retired
                            ! [GR-CROP C11] nmrain/rainamount/rainfluxarray/raintimearray retired from import:
                            !   readmeteo now writes directly to state%atmosphere%X
      ! [SS-GR-CROP A16] nutrient legacy globals — in WSN modules, not variables.f90
      use Wofost_Soil_Declarations, only: FOM_t, Bio_t, Hum_t, FOM_t0, Bio_t0, Hum_t0, &
                                          cNH4_t, cNO3_t, cNH4_t0, cNO3_t0, cNH4_av, cNO3_av, &
                                          Nminer, Cdissi, NsupplyNH4N, NsupplyNO3N, &
                                          FOM_old, Bio_old, Hum_old, NFOM_old, NBio_old, NHum_old, &
                                          NH4_old, NO3_old, &
                                          FOM_end, Bio_end, Hum_end, NFOM_end, NBio_end, NHum_end, &
                                          NH4_end, NO3_end, &
                                          FOM_add, NFOM_add, Hum_add, NHum_add, &
                                          FOM_cres, NFOM_cres, Hum_cres, NHum_cres, &
                                          NH4_intop, NH4_inlat, NH4_inbot, NH4_upt, NH4_out, NH4_nitrif, NH4_miner, &
                                          NO3_intop, NO3_inlat, NO3_inbot, NO3_upt, NO3_out, NO3_denitr, &
                                          FOM2Bio, NFOM2Bio, FOM2Hum, NFOM2Hum, FOM_dis, NFOM_min, &
                                          Bio2Bio, NBio2Bio, Bio2Hum, NBio2Hum, Bio_dis, NBio_min, &
                                          Hum2Bio, NHum2Bio, Hum2Hum, NHum2Hum, Hum_dis, NHum_min, &
                                          NH4N_amend, NO3N_amend, NH4N_cres, NO3N_cres, NH4N_volat, &
                                          iNLOSSL_1, iNLOSSR_1, iNLOSSS_1, iNLOSSO_1, &
                                          idwrt_1, idwlv_1, idwst_1, idwso_1, &
                                          Ntotuptake, Ptotuptake, &
                                          DMcressur, Ncressurf, Pcressurf, &
                                          DMcresbott, Ncresbott, Pcresbott
      use Wofost_Soil_Interface,    only: NdemandSoil, NsupplySoil, Ndemand, Nsupply, LaiCritNupt
      use soilwater_state_mod, only: soilwater_init
      use tillage_state_mod, only: tillage_init
      use drainage_mod, only: drainage_init
      use surfacewater_mod, only: SurfaceWater
      use tillage_mod,   only : DoTillage
      use swap_log, only: log_info
      use runoff_mod, only: cn_init
      use temperature_mod, only: Temperature
      use solute_mod, only: solute, solute_init
      use agetracer_mod, only: AgeTracer
      use soilgrid_mod, only: CalcGrid
      use soilhydraulics_mod, only: soilwater
      use config_to_variables_mod, only: config_to_variables
      use irrigation_mod, only: SSDI_irrigation
      use timecontrol_mod, only: timecontrol_init, itertime_init
      type(swap_state_t),          intent(out)   :: state
      type(swap_config_t), target, intent(inout) :: config  ! target: crop_config_global => config%crop (Task 3); inout: already loaded by caller
      logical :: request_smaller_dt   ! intent(out) dummy for SurfaceWater(1)

!  Initialization of all variables in Module Variables
   call Initialize
   ! [SS-GR-CROPRT A5] mirror flCropHarvest zero-init (Initialize has no state arg)
   state%crop%common%flCropHarvest = flCropHarvest

   ! [GR-CROP-DVS] non-owning config pointer; lifetime matches state's.
   ! Top-level compute routines can now read switches via state%cfg%X%Y
   ! without needing a separate config arg.
   state%cfg => config

!  iteration and timing statistics
   call itertime_init(state)

!  config_to_variables seeds state%timecontrol from config (Task 5).
   call config_to_variables(config, state)

   ! [GR-FINAL C1] state%timecontrol%iyear/imonth/dt now written directly by
   ! config_to_variables (tc_*_init_buf buffers retired).

!  shared simulation
   if (flSwapShared) call SharedSimulation(1)

!  initialize time variables and switches/flags
   call timecontrol_init(state)

!  calculate grid parameters
   ! [GR-BH Task 35] CalcGrid now writes directly to state%mesh%X; state%mesh%init bridge retired
   call CalcGrid(state)
   call soilwater_init(state%soilwater, state%mesh%numnod, numlay)   ! SS-CRP Phase 1 C-1.2: allocate per-node arrays + mfluxtable
   call state%nutrients%init(numlay)                                ! [SS-GR-CROP A11] zero nutrients state

   ! [SS-GR-CROP A16] dual-write nutrients — WSN organic matter pools + N coupling
   ! Primary state pools (scalar or maxfn=8 array)
   state%nutrients%fom_t   = FOM_t
   state%nutrients%bio_t   = Bio_t
   state%nutrients%hum_t   = Hum_t
   state%nutrients%fom_t0  = FOM_t0
   state%nutrients%bio_t0  = Bio_t0
   state%nutrients%hum_t0  = Hum_t0
   ! Mineral N concentrations
   state%nutrients%cnh4_t  = cNH4_t
   state%nutrients%cnh4_t0 = cNH4_t0
   state%nutrients%cnh4_av = cNH4_av
   state%nutrients%cno3_t  = cNO3_t
   state%nutrients%cno3_t0 = cNO3_t0
   state%nutrients%cno3_av = cNO3_av
   ! Per-timestep rates / scratchpads
   state%nutrients%nminer      = Nminer
   state%nutrients%cdissi      = Cdissi
   state%nutrients%nsupplynh4n = NsupplyNH4N
   state%nutrients%nsupplyno3n = NsupplyNO3N
   ! Crop–soil N coupling (from Wofost_Soil_Interface)
   state%nutrients%ndemandsoil  = NdemandSoil
   state%nutrients%nsupplysoil  = NsupplySoil
   state%nutrients%ndemand      = Ndemand
   state%nutrients%nsupply      = Nsupply
   state%nutrients%laicritnupt  = LaiCritNupt
   ! Balance-check state: OM
   state%nutrients%fom_old  = FOM_old
   state%nutrients%bio_old  = Bio_old
   state%nutrients%hum_old  = Hum_old
   state%nutrients%fom_end  = FOM_end
   state%nutrients%bio_end  = Bio_end
   state%nutrients%hum_end  = Hum_end
   state%nutrients%fom_add  = FOM_add
   state%nutrients%fom_cres = FOM_cres
   state%nutrients%fom2bio  = FOM2Bio
   state%nutrients%fom2hum  = FOM2Hum
   state%nutrients%fom_dis  = FOM_dis
   state%nutrients%bio2bio  = Bio2Bio
   state%nutrients%bio2hum  = Bio2Hum
   state%nutrients%bio_dis  = Bio_dis
   state%nutrients%hum_add  = Hum_add
   state%nutrients%hum_cres = Hum_cres
   state%nutrients%hum2bio  = Hum2Bio
   state%nutrients%hum2hum  = Hum2Hum
   state%nutrients%hum_dis  = Hum_dis
   ! Balance-check state: organic N
   state%nutrients%nfom_old  = NFOM_old
   state%nutrients%nbio_old  = NBio_old
   state%nutrients%nhum_old  = NHum_old
   state%nutrients%nfom_end  = NFOM_end
   state%nutrients%nbio_end  = NBio_end
   state%nutrients%nhum_end  = NHum_end
   state%nutrients%nfom_add  = NFOM_add
   state%nutrients%nfom_cres = NFOM_cres
   state%nutrients%nfom2bio  = NFOM2Bio
   state%nutrients%nfom2hum  = NFOM2Hum
   state%nutrients%nfom_min  = NFOM_min
   state%nutrients%nbio2bio  = NBio2Bio
   state%nutrients%nbio2hum  = NBio2Hum
   state%nutrients%nbio_min  = NBio_min
   state%nutrients%nhum_add  = NHum_add
   state%nutrients%nhum_cres = NHum_cres
   state%nutrients%nhum2bio  = NHum2Bio
   state%nutrients%nhum2hum  = NHum2Hum
   state%nutrients%nhum_min  = NHum_min
   ! Balance-check state: NH4
   state%nutrients%nh4_old    = NH4_old
   state%nutrients%nh4_end    = NH4_end
   state%nutrients%nh4_miner  = NH4_miner
   state%nutrients%nh4_intop  = NH4_intop
   state%nutrients%nh4_inlat  = NH4_inlat
   state%nutrients%nh4_inbot  = NH4_inbot
   state%nutrients%nh4_upt    = NH4_upt
   state%nutrients%nh4_out    = NH4_out
   state%nutrients%nh4_nitrif = NH4_nitrif
   ! Balance-check state: NO3
   state%nutrients%no3_old    = NO3_old
   state%nutrients%no3_end    = NO3_end
   state%nutrients%no3_intop  = NO3_intop
   state%nutrients%no3_inlat  = NO3_inlat
   state%nutrients%no3_inbot  = NO3_inbot
   state%nutrients%no3_upt    = NO3_upt
   state%nutrients%no3_out    = NO3_out
   state%nutrients%no3_denitr = NO3_denitr
   ! Amendment / residue N tracking
   state%nutrients%nh4n_amend = NH4N_amend
   state%nutrients%no3n_amend = NO3N_amend
   state%nutrients%nh4n_cres  = NH4N_cres
   state%nutrients%no3n_cres  = NO3N_cres
   state%nutrients%nh4n_volat = NH4N_volat
   ! Previous-period crop residue DM and N losses
   state%nutrients%idwrt_1   = idwrt_1
   state%nutrients%idwlv_1   = idwlv_1
   state%nutrients%idwst_1   = idwst_1
   state%nutrients%idwso_1   = idwso_1
   state%nutrients%inlossl_1 = iNLOSSL_1
   state%nutrients%inlossr_1 = iNLOSSR_1
   state%nutrients%inlosss_1 = iNLOSSS_1
   state%nutrients%inlosso_1 = iNLOSSO_1
   ! ANIMO cropext accumulators
   state%nutrients%ntotuptake = Ntotuptake
   state%nutrients%ptotuptake = Ptotuptake
   state%nutrients%dmcressur  = DMcressur
   state%nutrients%ncressurf  = Ncressurf
   state%nutrients%pcressurf  = Pcressurf
   state%nutrients%dmcresbott = DMcresbott
   state%nutrients%ncresbott  = Ncresbott
   state%nutrients%pcresbott  = Pcresbott

   ! [SS-GR-BH A6] soilwater layer flats — placed here because soilwater_init
   ! allocates the state arrays (nlay-sized) AFTER config_to_variables runs.
   ! All layer flats sourced directly from config (legacy globals retired Task 36).
   if (allocated(config%soil%hydraulics%ksatexm)) &
      state%soilwater%ksatexm(:) = config%soil%hydraulics%ksatexm(1:size(state%soilwater%ksatexm))
   if (allocated(config%soil%hydraulics%ksatfit)) &
      state%soilwater%ksatfit(:) = config%soil%hydraulics%ksatfit(1:size(state%soilwater%ksatfit))
   ! [GR-BH Task 36] cofani multi-source: drain.cofani first, soil.cofani overrides (soil wins).
   if (allocated(config%drain%cofani)) &
      state%soilwater%cofani(1:size(config%drain%cofani)) = config%drain%cofani
   if (allocated(config%soil%cofani)) &
      state%soilwater%cofani(1:size(config%soil%cofani))  = config%soil%cofani
   state%soilwater%flksatexm   = .false.   ! never set in adapter
   ! [GR-BH Task 36] orgmat multi-source: soil.orgmat first; heat.porg backfills when absent.
   if (allocated(config%soil%orgmat)) &
      state%soilwater%orgmat(1:size(config%soil%orgmat))  = config%soil%orgmat
   if (.not. allocated(config%soil%orgmat) .and. allocated(config%heat%porg)) &
      state%soilwater%orgmat(1:min(size(config%heat%porg), size(state%soilwater%orgmat))) = &
         config%heat%porg(1:min(size(config%heat%porg), size(state%soilwater%orgmat)))
   if (allocated(config%heat%psand)) &
      state%soilwater%psand(:) = config%heat%psand(1:size(state%soilwater%psand))
   if (allocated(config%heat%psilt)) &
      state%soilwater%psilt(:) = config%heat%psilt(1:size(state%soilwater%psilt))
   if (allocated(config%heat%pclay)) &
      state%soilwater%pclay(:) = config%heat%pclay(1:size(state%soilwater%pclay))
   ! [SS-GR-BH A7] seed soilwater runtime scalars; swbotb_runtime sourced from config.
   state%soilwater%swbotb_runtime = config%bottom_boundary%swbotb
   state%soilwater%q0    = 0.0d0
   state%soilwater%k1max = 0.0d0
   state%soilwater%H0max = 0.0d0
   ! [GR-FINAL C1/C2] seed state%soilwater/atmosphere/timecontrol from config.
   ! swinco=3 warm-restart: pond/pondini/dt/h-profile/atmosphere from soil.initial.
   ! swinco<3: pondini/pond from soil.pondini; atmosphere inits to zero.
   ! Both blocks consolidated here to eliminate the duplicate swinco==3 guard.
   call state%atmosphere%init(config)                     ! GR-ATM: zero flat scalars + cohorts + snapshot config-derived params (snowcoef/swsublim/swetsine)
   if (config%soil%swinco == 3 .and. &
       allocated(config%soil%initial%h_file) .and. &
       len_trim(config%soil%initial%h_file) > 0) then
      ! soilwater init from warm-restart record
      state%soilwater%pondini = config%soil%initial%pond
      state%soilwater%pond    = config%soil%initial%pond
      ! dt from warm-restart (soil.initial.dt supersedes simulation%numerical%dt for swinco=3)
      state%timecontrol%dt = config%soil%initial%dt
      ! h profile: CSV re-read after soilwater_init has allocated state%soilwater%h
      block
         use csv_reader_mod, only: read_csv_table
         use error_mod,      only: error_collection_t
         real(8), allocatable     :: tbl(:,:)
         type(error_collection_t) :: errs
         character(len=2)         :: hdr(2)
         integer :: nrows, ki
         hdr(1) = 'z '
         hdr(2) = 'h '
         call read_csv_table(trim(config%soil%initial%h_file), hdr, tbl, errs)
         call errs%abort_if_fatal()
         nrows = size(tbl, 1)
         do ki = 1, min(nrows, size(state%soilwater%h))
            state%soilwater%h(ki) = tbl(ki, 2)
         end do
      end block
      ! atmosphere warm-restart [SS-ATM A-2.6]: ssnow/ldwet/slw from soil.initial
      ! spev/saev not in config; remain zero from atmosphere%init
      state%atmosphere%ssnow = config%soil%initial%ssnow
      state%atmosphere%ldwet = config%soil%initial%ldwet
      state%atmosphere%slw   = config%soil%initial%slw
      if (config%meteo%snow%swsnow /= 1) state%atmosphere%ssnow = 0.0d0
   else
      state%soilwater%pondini = config%soil%pondini
      state%soilwater%pond    = config%soil%pondini   ! legacy alias: pond <-> pondini for swinco<3
   end if

   ! [GR-CROP-DVS] atmosphere meteo + runoff CN seeds retired:
   ! all legacy globals (arad/atmn/atmx/ahum/awin/arai/aetr/wet/atav/epot/tpot/grain/nrain/
   ! tav/daynrfirst/daynrlast/atmin7/nofd/isua/avevaptb/avprectb/pfreetb/pstemtb/scanopytb/
   ! CNdry/CNwet/ThetaRef/Runoff_CN/wc_cor/wc10/iCNtab/CNtimTAB/CNrefTAB) were zero at init;
   ! state fields default-init to zero; readmeteo writes state directly at runtime.

   ! [GR-CROP C11] rain timing dual-write retired: readmeteo now writes directly to state%atmosphere%X

   ! [GR-ATM C8] out_tmn/tmx/hum/win/etr/wet/rad seeding dropped: legacy globals retired; state%atmosphere%X written by meteoday

   ! [SS-GR-ATM A12] seed state%crop from legacy crop globals; lai retired
   state%crop%swcfbs          = swcfbs
   state%crop%flCropEmergence = flCropEmergence
   ! [SS-GR-CROP A14] dual-write crop_common; dvs/tsum retired (writes target state directly)
   state%crop%common%daycrop        = daycrop
   state%crop%common%swcrp          = swcrp
   state%crop%common%icrop          = icrop
   state%crop%common%flCropCalendar = flCropCalendar
   state%crop%common%flCropOutput   = flCropOutput
   state%crop%common%flCropNut      = flCropNut
   state%crop%common%flHarvestDay   = flHarvestDay
   ! [SS-GR-CROPRT C1] swend dual-write dropped — state%crop%common%swend retired; ADR 0009: always 0
   ! [SS-GR-CROP A15] dual-write crop_wofost
   state%crop%wofost%swbulb   = (swbulb == 1)   ! integer→logical conversion
   ! [SS-GR-CROP A15] dual-write crop_grass
   ! seqgrazmow/seqgrazmowpot/dateharvest dual-writes retired — see state%crop%grass
   state%crop%grass%swpotrelmf    = swpotrelmf
   ! [SS-GR-CROP A15] dual-write crop_fixed retired — see state%crop%fixed%X

   ! [SS-TC TC-14] alias TC fields used in init block
   block
   associate( &
      flSnow         => state%timecontrol%flSnow,         &
      flSolute       => state%timecontrol%flSolute,       &
      flSurfaceWater => state%timecontrol%flSurfaceWater, &
      flTemperature  => state%timecontrol%flTemperature )

   call tillage_init(state%tillage, numlay)         ! SS-TIL T-2: allocate/zero tillage state unconditionally
   if (flTillage) call DoTillage(1, state)
   if (flSSDI)    call SSDI_irrigation(1, state)  ! [SS-SWC S-2.12B]

!  Allocate and initialise heat state arrays before SoilWater(1) so that
!  hconduc can read state%heat%tsoil(node) during hydraulic-conductivity init.
!  SS-SWC S-2.2: heat_init moved earlier to satisfy mandatory tsoil_node arg.
   call state%heat%init(config%heat, state%mesh%numnod)  ! GR-BH Task 11: type-bound init

!  initialize SoilWater rate/state variables
   call SoilWater(1, state)
   ! SS-ATM A-2.6: state added — CNmethod signature updated for retired nraidt/melt
   if (swuseCN == 1) call cn_init(state)

!  Allocate and initialise drainage state arrays.  Config is passed so
!  drainage_init can seed state%drainage%wetper(1) from config%drain%wetper
!  (dramet==2) without reading the now-deleted legacy global wetper.
!  ADR 0031 Phase 2 Task 5: drainl/wetper/ztopdislay/qdrd globals deleted.
   call drainage_init(state, config)
   ! [SS-GR-BH A9 / GR-BH Task 37] state%drainage scalars + geometry sourced directly from config.
   ! Bare globals nrlevs/swdivd/swnrsrf/swtopnrsrf/swdivdinf/FacDpthInf/L/zbotdr deleted from variables.f90.
   state%drainage%nrlevs     = config%drain%nrlevs
   state%drainage%swdivd     = config%drain%swdivd
   state%drainage%swnrsrf    = config%drain%surface_runoff%swnrsrf
   state%drainage%swtopnrsrf = config%drain%surface_runoff%swtopnrsrf
   state%drainage%swdivdinf  = config%drain%surface_runoff%swdivdinf
   state%drainage%FacDpthInf = config%drain%surface_runoff%facdpthinf
   ! L: multi-source — dramet==2 uses scalar lm (m→cm); per-level uses config%drain%L(:) (already in cm).
   if (config%drain%dramet == 2) then
      state%drainage%L(1) = 100.0d0 * config%drain%lm
   else if (allocated(config%drain%L)) then
      state%drainage%L(1:size(config%drain%L)) = config%drain%L
   end if
   ! zbotdr: multi-source — dramet==2 uses scalar zbotdr_basic; per-level uses config%drain%zbotdr(:).
   if (config%drain%dramet == 2) then
      state%drainage%zbotdr(1) = config%drain%zbotdr_basic
   else if (allocated(config%drain%zbotdr)) then
      state%drainage%zbotdr(1:size(config%drain%zbotdr)) = config%drain%zbotdr
   end if
   ! owltab: populated via CSV loop in config_to_variables.f90 — still uses bare global.
   ! Deferred: blocked by nowltab which is not yet in state (runtime reads nowltab from variables module).
   state%drainage%owltab(:,:) = owltab(1:size(state%drainage%owltab,1), &
                                       1:size(state%drainage%owltab,2))
   if (flSolute) call solute_init(state)   ! SS-SLST Phase 2 Task 7: seed state%solute from config-populated globals

!  initialize SurfaceWater management variables
   ! S4 state init for surfacewater (spec docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).
   if (flSurfaceWater) call state%surfacewater%init(config%surface_water, config%drain, state%mesh%numnod)
   ! SurfaceWater(task=1)'s case(1) is now a no-op stub; init was hoisted to the line above.
   ! Dispatcher-case removal is a separate cleanup follow-up.
   if (flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)

!  [MACRO-RETIRE 2026-05-12] MACROPORE call retired (ADR 0040). flMacroPore is permanent .false.
!  if (flMacroPore) call MACROPORE(1, state)

!  initialize SoilTemperature rate/state variables
   if (flTemperature) call Temperature(1, state, config)

!  initialize Snow: handshake snowinco<->ssnow based on swinco (formerly snow_init)
   if (flSnow) then
      if (config%soil%swinco == 3) then
         state%atmosphere%snowinco = state%atmosphere%ssnow
      else
         state%atmosphere%ssnow = state%atmosphere%snowinco
      end if
   end if

!  initialize Solute rate/state variables
   if (flSolute) call Solute(1, state)

!  [SS-GR-CROPRT A1] AgeTracer init dropped — flAgeTracer retired (ADR 0032; body always unreachable)

!  Soil Management init: SoilManagement(1) was the legacy reader entry
!  point and is now a no-op (SS-C step 3). flCropNut is now driven by
!  the per-rotation typed config (ADR 0028); the SoilManagement(2..7)
!  call sites below run when a rotation has flcropnut=true.

!  open Output files and write headers (always run; iCaller branch retired)
   call SwapOutput(1, state)
   call SoilWaterOutput(1, state, config)   ! [SS-GR-CROPRT A3] config threaded for swcsv/swcsv_tz
!  ADR 0009 Phase 5+: IrrigationOutput deleted (swirg=0).
   if (flTemperature)  call TemperatureOutput(1, state)
   if (flSolute)       call SoluteOutput(1, state)
   ! [SS-GR-CROPRT A1] AgeTracerOutput(1) dropped — flAgeTracer retired (ADR 0032)
   if (flSnow)         call SnowOutput(1, state)
   ! [MACRO-RETIRE 2026-05-12] MacroPoreOutput retired (ADR 0040).
   if (flSurfaceWater) call SurfaceWaterOutput(1, state)
   end associate
   end block

   call log_info('swap', 'Initialization complete for project: ' // trim(project))

   end subroutine swap_init_from_loaded_config

   subroutine swap_run_step(state, config)
      use variables, only : flswapshared, flcropnut, swfrost, &
                            ! [SS-GR-CROPRT A1] flagetracer dropped — retired (ADR 0032; always .false.)
                            flcropcalendar, &
                            flharvestday, flcropoutput, swcrp, &
                            ! [GR-FINAL C3] swend dropped: read via state%crop%common%swend (C category, inv. §C)
                            flTillage, flSSDI
      use cropgrowth_helpers_mod, only: CropOutput  ! GR-CROPWS Phase 0
      use timestep_control_mod, only: fldecdt
      use timecontrol_mod, only: timecontrol_advance, timecontrol_reduce_dt, &
                                  timecontrol_day_end, itertime_check
      use surfacewater_mod, only: SurfaceWater, surfacewater_year_reset
      use tillage_mod, only: DoTillage
      use boundbottom_mod, only: BoundBottom
      use meteo_mod, only: ProcessMeteoDay
      use meteo_process_mod, only: ReadMeteoDay
      use snow_mod, only: snow_step
      use meteodt_mod, only: MeteoDT
      use rootextraction_mod, only: RootExtraction
      use frozencond_mod, only: FrozenCond, FrozenBounds
      use temperature_mod, only: Temperature
      use solute_mod, only: solute
      use agetracer_mod, only: AgeTracer
      use soilhydraulics_mod, only: soilwater, SoilWaterStateVar
      use irrigation_mod, only: irrigation, SSDI_irrigation
      use management_soil_mod, only: SoilManagement
      use drainage_mod, only: drainage
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      logical :: request_smaller_dt
      logical, external :: dtleap

      interface
         subroutine CropGrowth(task, tsoil, state)
            use swap_state_mod, only: swap_state_t
            integer, intent(in) :: task
            real(8), intent(in) :: tsoil(:)
            type(swap_state_t), intent(inout) :: state
         end subroutine CropGrowth
      end interface

!  [SS-TC TC-14] bind TC aliases for all timestep-loop fields used below
   associate( &
      tc_flYearStart => state%timecontrol%flYearStart, &
      tc_flDayStart  => state%timecontrol%flDayStart,  &
      tc_flDayEnd    => state%timecontrol%flDayEnd,    &
      tc_daynr       => state%timecontrol%daynr,       &
      tc_iyear       => state%timecontrol%iyear,       &
      flrunend       => state%timecontrol%flRunEnd,    &
      flOutput       => state%timecontrol%floutput,    &
      flOutputShort  => state%timecontrol%floutputshort,&
      flMeteoDt      => state%timecontrol%flmeteodt,   &
      flETSine       => state%timecontrol%fletsine,    &
      flSnow         => state%timecontrol%flSnow,      &
      flSolute       => state%timecontrol%flSolute,    &
      flTemperature  => state%timecontrol%flTemperature,&
      flDrain        => state%timecontrol%flDrain,     &
      flSurfaceWater => state%timecontrol%flSurfaceWater,&
      flIrrigate     => state%timecontrol%flIrrigate,  &
      fldtreduce     => state%timecontrol%fldtreduce )

!     get Meteo data
   if (tc_flYearStart) call ReadMeteoYear(state, config)  ! SS-TC TC-13 / SS-GR-FINAL B1

      if (tc_flDayStart) then  ! SS-TC TC-13

!        read meteo data for current day
         call ReadMeteoDay(state, config)  ! SS-ATM A-1.6 / SS-GR-ATM B23: state + config threaded

!        check growing season
         call CropGrowth(1, state%heat%tsoil, state)  ! SS-CRP C-1.3: state added for dual-write

!        calculate Irrigation rate/state variables
         if (flIrrigate) call irrigation(2, state)

!        process Meteo data
         call ProcessMeteoDay(state, config)  ! SS-GR-ATM B24: config added
         if (flTillage) call DoTillage(2, state)

      end if

!     process Meteo data
      if (flMeteoDt .or. flETSine) call MeteoDT(state)

!     shared simulation
      if (flSwapShared .and. tc_flDayStart) call SharedSimulation(2)  ! SS-TC TC-13

!     calculate Snow: MH+MM - probably to be moved within IF-block above, prior to call ProcessMeteoDay ...
      ! SS-HEAT Phase 2 Task 6: pass state so Snow reads tsoil from state%heat
      if (flSnow .and. tc_flDayStart) call snow_step(state)  ! SS-TC TC-13

!     calculate reduction for conductivities for frozen conditions
      if (SwFrost.eq.1) then
         call FrozenCond(state, config)
      end if

!     calculate potential and actual root water extraction profile
      call RootExtraction(state)

!     determine SoilWater bottom boundary conditions
      call BoundBottom(state, config)  ! [SS-GR-BH Task 18]: config added for bottom_boundary fields

      fldtreduce = .true.
      do while(fldtreduce)
         fldtreduce = .false.

!        calculate drainage fluxes
         if (fldrain)                           call Drainage(state)
         ! SS-SWST Phase 2: SurfaceWater sets request_smaller_dt; propagate to fldecdt here.
         if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(2, state, request_smaller_dt)
         if (request_smaller_dt) fldecdt = .true.
         if (SwFrost.eq.1)                      call FrozenBounds(state, config)

!        calculate SoilWater, incl macropores (headcalc inside may also set fldecdt on non-convergence)
         if (.not.fldecdt) call SoilWater(2, state)

!        calculate surface water balance
         if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(3, state, request_smaller_dt)
         if (request_smaller_dt) fldecdt = .true.

!        update time variables and switches/flags
         if (fldecdt) then
            call SoilWaterStateVar(2, state)
            call timecontrol_reduce_dt(state)
            fldtreduce = .true.
         end if

      end do

!     calculate SoilWater rate/state variables
      call SoilWater(3, state)

!     calculate SoilTemperature rate/state variables
   if (flTemperature) call Temperature(2, state, config)

!     calculate Solute rate/state variables
      if (flSolute) call Solute(2, state)

!     [SS-GR-CROPRT A1] AgeTracer(2) dropped — flAgeTracer retired (ADR 0032)

!     update time variables and switches/flags
      call timecontrol_advance(state)

!     at the end of a day,
      if (tc_flDayEnd) then  ! SS-TC TC-13

!        update Soil nutrient status variables
         if (flCropNut) call SoilManagement(2, state)

!        calculate potential crop growth
         if (flCropCalendar) call CropGrowth(2, state%heat%tsoil, state)

!        amendent of crop residues from previous day
         if (flCropNut) call SoilManagement(5, state)

!        amendent of fertilizers of current day
         if (flCropNut) call SoilManagement(3, state)

!        calculate actual crop growth (calculation of actual crop rate and state variables)
         if (flCropCalendar) call CropGrowth(3, state%heat%tsoil, state)

!        Simulate Soil Nutrient processes
         if (flCropNut) call SoilManagement(4, state)

!        harvest of crop
         if (flCropCalendar) call CropGrowth(4, state%heat%tsoil, state)

!        timing statistics : prevent (near) endless simulations
         if (state%timecontrol%flMaxIterTime) call itertime_check(state)  ! [SS-BMI2 Task 4]

!        Better here: check if subsurface irrigation is required for next day,
!                     and determine if time step needs to be changed due to dt_SSDI_event
         if (flSSDI) call SSDI_irrigation(2, state)  ! [SS-SWC S-2.12B]
         call timecontrol_day_end(state)

      end if

!     output section
         if (flOutput) then
            call SwapOutput(2, state)
            call SoilWaterOutput(2, state, config)   ! [SS-GR-CROPRT A3]
            if (flTillage) call DoTillage(3, state)
            if (flTemperature)   call TemperatureOutput(2, state)
            if (flSolute)        call SoluteOutput(2, state)
            ! [SS-GR-CROPRT A1] AgeTracerOutput(2) dropped — flAgeTracer retired (ADR 0032)
            if (flSnow)          call SnowOutput(2, state)
            ! [MACRO-RETIRE 2026-05-12] MacroPoreOutput retired (ADR 0040).
            if (flSurfaceWater) then
               if (tc_daynr == merge(366, 365, dtleap(tc_iyear))) &  ! SS-TC TC-13
                  call surfacewater_year_reset(state%surfacewater)
               call SurfaceWaterOutput(2, state)
            end if
         else
            if (flOutputShort)   call SoilWaterOutput(2, state, config)   ! [SS-GR-CROPRT A3]
         end if
         if (tc_flDayEnd .and. (flOutput .or. flHarvestDay)) then  ! SS-TC TC-13
            if (flCropCalendar .and. flCropOutput) then
               if (swcrp.eq.1) call CropOutput(2, state)
            end if
         end if
!        ADR 0009 Phase 5+: IrrigationOutput deleted (swirg=0).
         if (tc_flDayEnd .and. flCropNut)    call SoilManagement(6, state)   ! SS-TC TC-13
         ! [SS-GR-CROPRT C1] swend.eq.2 daily-dump branch dropped — ADR 0009: swend always 0

!    shared simulation
     if (flSwapShared .and. tc_flDayEnd) call SharedSimulation(3)  ! SS-TC TC-13

   end associate  ! SS-TC TC-13: tc_flYearStart, tc_flDayStart, tc_flDayEnd, tc_daynr, tc_iyear

   end subroutine swap_run_step

   subroutine swap_close(state, config)
      use variables, only : flswapshared, flcropnut, project, swcrp
                            ! [SS-GR-CROPRT A1] flagetracer dropped — retired (ADR 0032; always .false.)
                            ! [SS-GR-CROPRT C1] swend dropped — global + state field retired (ADR 0009)
      use swap_log,  only: log_info
      use management_soil_mod, only: SoilManagement
      use timecontrol_mod, only: itertime_close
      use cropgrowth_helpers_mod, only: CropOutput  ! GR-CROPWS Phase 0
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config  ! unused: kept for parallel signature with swap_init/swap_run_step

!  iteration and timing statistics
   call itertime_close(state)

!  close output files (always run; iCaller branch retired)
   if (flSwapShared) call SharedSimulation(4)
   call SwapOutput(3, state)
   ! [SS-GR-CROPRT C1] swend.eq.1 end-sim-dump branch dropped — ADR 0009: swend always 0
   call SoilWaterOutput(4, state, config)   ! [SS-GR-CROPRT A3]
   if (swcrp.eq.1) call CropOutput(3, state)
   ! [SS-TC TC-14] flag reads via state%timecontrol
   if (state%timecontrol%flTemperature)  call TemperatureOutput(3, state)
   if (state%timecontrol%flSolute)       call SoluteOutput(3, state)
   ! [SS-GR-CROPRT A1] AgeTracerOutput(3) dropped — flAgeTracer retired (ADR 0032)
!  ADR 0009 Phase 5+: IrrigationOutput deleted (swirg=0).
   if (state%timecontrol%flSnow)         call SnowOutput(3, state)
   ! [MACRO-RETIRE 2026-05-12] MacroPoreOutput retired (ADR 0040).
   if (state%timecontrol%flSurfaceWater) call SurfaceWaterOutput(3, state)
   if (flCropNut)                        call SoilManagement(7, state)

!  write okay file for external use
   call WriteSwapOk(Project)

   call log_info('swap', 'Simulation complete for project: ' // trim(project))

   end subroutine swap_close

end module swap_mod
