!> @file legacy_state.f90
!! Sweep 3 transitional bag: holds globals that were in variables.f90.
!! Generated for retirement arc — single struct lets us delete variables.f90.
module legacy_state_mod
   use, intrinsic :: iso_fortran_env, only: real64, int32
   use swap_array_dimensions  ! pulls all MA* constants used by remaining decls
   implicit none
   private
   public :: legacy_state_t

   type :: legacy_state_t
      ! [GR-ATM 2026-05-23] logf retired — swap_log owns the log-file unit
      ! [GR-CROP 2026-05-25] o2_* workspace cluster retired — all 19 fields
      ! migrated to state%crop%oxygen (w_root, w_root_z0, soil_temp,
      ! sat_water_cont, gas_filled_porosity, d_o2inwater, d_root, d_soil,
      ! perc_org_mat, soil_density, depth, shape_factor_microbialr,
      ! root_radius, r_microbial_z0, waterfilm_thickness, bunsencoeff,
      ! c_min_micro, c_macro, ctopnode).
      real(real64), allocatable :: outdat(:)  !! Array with output dates for water and solute balances
      real(real64), allocatable :: outdatint(:)  !! Array with intermediate output dates
      logical :: flSwapShared  !! Flag to indicate the shared simultaneous simulation with other applications
      character(len=16) :: outfil  !! Name of output file
      character(len=80) :: pathwork  !! Path to work directory
      character(len=80) :: project  !! Name of project
      integer :: daynrfirst  !! First calendar day number for which meteorological data is available in current year
      integer :: daynrlast  !! Last calendar day number for which meteorological data is available in current year
      ! [GR-ATM 2026-05-23] detrecord/irectotal retired — see state%atmosphere%{detrecord,irectotal}
      ! [GR-CROP 2026-05-25] nofd retired — see state%atmosphere%nofd
      integer :: swetsine  !! Switch: 0 = Tp and Ep uniform during a day; 1 = Tp and Ep are distributed as sine waves during a day
      integer, allocatable :: ad(:)  !! Array with day numbers in meteo file
      integer, allocatable :: am(:)  !! Array with month numbers in meteo file
      real(real64), allocatable :: aetr(:)  !! Array with daily ETref input data (L/T)
      real(real64), allocatable :: ahum(:)  !! Array with daily humidity input data (M/L/T2)
      ! [GR-ATM 2026-05-23] angstroma/b retired — see state%atmosphere%angstrom{a,b}
      real(real64), allocatable :: arad(:)  !! Array with daily radiation input data (M/T2)
      real(real64), allocatable :: arai(:)  !! Array with daily precipitation sum input data (L/T)
      real(real64), allocatable :: atav(:)  !! In case of detailed weather input, air temperature of each weather record (L/T)
      real(real64), allocatable :: atmin7(:)  !! Array with minimum temperatures of last week (C)
      real(real64), allocatable :: atmn(:)  !! Array with daily minimum temperature input data (  )
      real(real64), allocatable :: atmx(:)  !! Array with daily maximum temperature input data (  )
      real(real64), allocatable :: awin(:)  !! Array with daily wind speed input data (L/T)
      ! [GR-ATM 2026-05-23] daylp retired — see state%atmosphere%daylp
      ! [GR-ATM 2026-05-23] dethum/detrad/detrain/dettav/dettime/detwind retired — see state%atmosphere%det*
      ! [GR-ATM 2026-05-23] dtEventRain retired — see state%atmosphere%dtEventRain
      real(real64), allocatable :: epot(:)  !! In case of detailed weather input, calculated Epot of each weather record (L/T)
      ! [GR-ATM 2026-05-23] cfevappond retired — see state%atmosphere%cfevappond
      ! [GR-ATM 2026-05-23] finterception retired — see state%atmosphere%finterception
      real(real64), allocatable :: grain(:)  !! In case of detailed weather input, gross rain flux of each weather record (L/T)
      real(real64), allocatable :: nrain(:)  !! In case of detailed weather input, calculated netto rain of each weather record (L/T)
      ! [GR-ATM 2026-05-23] rad/tav/tmn retired — see state%atmosphere%{rad,Tav,tmn}
      real(real64) :: tmnr  !! Average of minimum air temperature during past 7 days (oC)
      ! [GR-ATM 2026-05-23] tmx retired — see state%atmosphere%tmx
      real(real64), allocatable :: tpot(:)  !! In case of detailed weather input, calculated Tpot of each weather record (L/T)
      real(real64), allocatable :: wet(:)  !! Fraction of each day the crop is wet (L)
      character(len=200) :: metfil  !! Name of meteorological input file
      character(len=80) :: pathatm  !! Path to folder with meteorological input files
      ! [GR-ATM 2026-05-23] tsunrise_atm/tsunset_atm retired — see state%atmosphere%tsun{rise,set}_atm
      integer :: nod10_cn  !! Node at -10cm for CN runoff method - from meteoday.f90 CNmethod
      integer :: icn_atm  !! Current position in CN time table - from meteoday.f90 CNmethod
      ! [GR-CROP 2026-05-25] irrigevent retired — runtime-local in src/crop/irrigation.f90
      ! [GR-CROP 2026-05-25] irtype retired — canonical home is state%crop%irrigation%irtype
      ! [GR-CROP 2026-05-25] isua retired — canonical home is state%atmosphere%isua
      ! [GR-CROP 2026-05-25] isuas retired — schedule==1 dead branch
      integer :: nirri  !! retained — still written by src/core/timecontrol_mod.f90 (initial reset to 1)
      integer :: swirfix  !! retained — still consumed by src/core/timecontrol_mod.f90 (sets flIrrigate)
      ! [GR-CROP 2026-05-25] swcirrthres retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] cirrs retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] cirrthres retired — schedule==1 dead branch
      real(real64) :: dcrit  !! Depth (L) of sensor for soil water pressure head or water content
      ! [GR-CROP 2026-05-25] ditab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] dwatab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] fidtab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] hcritab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] irconc retired — canonical home is state%crop%irrigation%irconc
      ! [GR-CROP 2026-05-25] irdate retired — canonical home is state%crop%irrigation%irdate
      ! [GR-CROP 2026-05-25] irdepth retired — canonical home is state%crop%irrigation%irdepth
      ! [GR-ATM 2026-05-23] nird retired — see state%atmosphere%nird
      ! [GR-CROP 2026-05-25] perirrsurp retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] raithreshold retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] rawtab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] tawtab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] tcritab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] tstairrig retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] tendirrig retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] treltab retired — schedule==1 dead branch
      ! [GR-CROP 2026-05-25] dayfix retired — schedule==1 dead branch (local-only counter)
      ! [GR-CROP 2026-05-25] tsumtime/tsumtemp/tsumdepth retired — orphan
      !   (only consumer was cropgrass_init use-list, never read in body)
      logical :: flCropCalendar  !! Flag indicating that crop season is active (but currently might be bare or cropped)
      logical :: flCropEmergence  !! Flag indicating period from crop emergence until harvest
      logical :: flCropHarvest  !! Flag indicating period from crop harvest until the end of crop season
      ! [GR-CROP 2026-05-25] flCropReadFile retired — see state%crop%common%flCropReadFile
      ! [GR-CROP 2026-05-25] flCropOpenFile retired — see state%crop%common%flCropOpenFile
      logical :: flCropOutput  !! Flag indicating writing of output.crp
      ! [GR-ATM 2026-05-23] croptype retired — see state%crop%common%croptype
      integer :: swcrp  !! Switch for output file *.CRP with daily crop output: 0 = no; 1 = yes
      integer :: crp  !! Internal number of crop output file *.CRP
      integer :: daycrop  !! Number of days that a crop exists
      integer :: icrop  !! Current crop number
      ! [GR-CROP 2026-05-25] idsl retired — see state%crop%wofost%idsl
      ! [GR-CROP 2026-05-25] noddrz retired — see crop_common_state_t%noddrz
      ! [GR-ATM 2026-05-23] atmtr retired — see state%atmosphere%atmtr
      ! [GR-CROP 2026-05-25] agerm/bgerm/cgerm retired — read direct from
      ! crop_config_global%rotation_wofost(icrop)%germination (cgerm/bgerm derived locally).
      ! [GR-CROP 2026-05-25] swpotrelmf retired — see state%crop%grass%swpotrelmf
      real(real64), allocatable :: avevaptb(:)  !! Gash interception model: average evaporation intensity during shower (-) as function of time (T)
      real(real64), allocatable :: avprectb(:)  !! Gash interception model: average rainfall intensity (-) as function of time (T)
      ! [GR-CROP 2026-05-25] c_mroot retired — see state%crop%oxygen%c_mroot
      real(real64), allocatable :: cropend(:)  !! Array with crop end dates
      real(real64), allocatable :: cropstart(:)  !! Array with crop start dates
      ! [GR-ATM 2026-05-23] difpp retired — see state%atmosphere%difpp
      ! [GR-CROP 2026-05-25] dlc retired — see state%crop%wofost%dlc
      ! [GR-CROP 2026-05-25] dlo retired — see state%crop%wofost%dlo
      ! [GR-ATM 2026-05-23] dsinbe retired — see state%atmosphere%dsinbe
      ! [GR-CROP 2026-05-25] f_senes retired — see state%crop%oxygen%f_senes
      ! [GR-CROP 2026-05-25] gasst/gasstpot retired — local SAVE in cropwofost_runtime%wofost
      ! [GR-ATM 2026-05-23] gc retired — see state%crop%common%gc
      ! [GR-CROP 2026-05-25] hdrygerm/hwetgerm retired — see cropwofost_config_t%germination.
      ! [GR-CROP 2026-05-25] max_resp_factor retired — see state%crop%oxygen%max_resp_factor
      ! [GR-CROP 2026-05-25] mrest/mrestpot retired — local SAVE in cropwofost_runtime%wofost
      ! [GR-CROP 2026-05-25] mrftb retired — always-zero on TOML (legacy reader removed)
      real(real64), allocatable :: pfreetb(:)  !! Gash interception model: free throughfall coefficient (-) as function of time (T)
      real(real64), allocatable :: pstemtb(:)  !! Gash interception model: stem flow coefficient (-) as function of time (T)
      ! [GR-ATM 2026-05-23] siccaptb retired — see state%crop%common%siccaptb
      ! [GR-CROP 2026-05-25] siccaplai retired — Gash interception (swinter=3) stub-errored on TOML
      ! [GR-CROP 2026-05-25] q10_microbial retired — see state%crop%oxygen%q10_microbial
      ! [GR-CROP 2026-05-25] q10_root retired — see state%crop%oxygen%q10_root
      ! [GR-CROP 2026-05-25] reltr retired — see state%crop%common%reltr
      ! [GR-CROP 2026-05-25] rid retired — see state%crop%common%rid
      ! [GR-CROP 2026-05-25] rootcoefa/rooteff/rootradius retired — JvL-only (swdrought=2);
      ! body extracted to src/crop/dormant/jongvanlier.f90.
      real(real64), allocatable :: scanopytb(:)  !! Gash interception model: storage capacity of canopy (-) as function of time (T)
      ! [GR-CROP 2026-05-25] shape_factor_rootr retired — see state%crop%oxygen%shape_factor_rootr
      ! [GR-CROP 2026-05-25] specific_resp_humus retired — see state%crop%oxygen%specific_resp_humus
      ! [GR-CROP 2026-05-25] tadw/tadwpot retired — local SAVE in cropwofost_runtime%wofost
      ! [GR-CROP 2026-05-25] tsumemeopt/TBASEM/TEFFMX retired — see cropwofost_config_t%germination.
      ! [GR-CROP 2026-05-25] tsumgerm retired — see crop_common_state_t%tsumgerm.
      ! [GR-CROP 2026-05-25] w_root_ss retired — see state%crop%oxygen%w_root_ss
      ! [GR-CROP 2026-05-25] wiltpoint retired — see state%crop%common%hlim4
      ! [GR-CROP 2026-05-25] wrtb retired — always-zero on TOML (legacy reader removed)
      ! [GR-CROP 2026-05-25] wrtmin retired — see state%crop%wofost%wrtmin
      ! [GR-CROP 2026-05-25] gwrt retired — see state%crop%wofost%gwrt
      ! [GR-CROP 2026-05-25] flanthesis retired — see state%crop%wofost%flanthesis
      logical :: flHarvestDay  !! Flag indicating that current day is harvest day
      character(len=40), allocatable :: cropfil(:)  !! Array with names of crop files
      ! [GR-CROP 2026-05-25] pathcrop retired — config%general%pathcrop is the canonical read.
      real(real64) :: rdmax  !! Maximum rooting depth in soil profile (L)
      ! [GR-ATM 2026-05-23] flCO2 retired — see state%atmosphere%flco2
      ! [GR-CROP 2026-05-25] co2amaxtb/co2efftb/co2tratb/co2year/co2ppm retired —
      ! see cropwofost_config%co2; flco2 dormant (no live consumer).
      ! [GR-CROP 2026-05-25] verndvs/vernsat/vernbase/vernrtb retired — see cropwofost_init_mod cw_vern*
      integer :: swbulb  !! switch to enable simulation of bulb crops (-)
      ! [GR-CROP 2026-05-25] drbl/drblpot/fbl retired — local SAVE in cropwofost_runtime%wofost
      ! [GR-CROP 2026-05-25] pld retired — see crop_config_global%rotation_wofost(icrop)%bulb%pld
      ! [GR-CROP 2026-05-25] remoc retired — see crop_config_global%rotation_wofost(icrop)%bulb%remoc
      ! [GR-CROP-DVS] plwt retired — see state%crop%wofost%plwt
      ! [GR-CROP 2026-05-25] flCropNut retired — see state%crop%common%flCropNut
      ! [GR-CROP 2026-05-25] nmxlv/ilnmxl/fstr retired — see cropwofost_init_mod cw_*
      ! [SS-GR-FINAL D1] amFERT retired — 0 consumers
      ! [GR-CROP 2026-05-25] till_* Group AB stubs retired — migrated to
      !   state%tillage (Ntill/Ntypes/i_n_model/iRedist/Max_Z_tillage/...).
      !   swtill canonical home is state%cfg%soil%swtill.
      ! [GR-SOIL 2026-05-24] iHWCKmodel retired — orphan stub; state%soilwater%iHWCKmodel.
      logical, allocatable :: BiModal(:)  !! logical indicating whether chosen model is bi-modal or not
      logical, allocatable :: NoVap(:)  !! logical indicating that NO vapour flow is to be considered in PDI K-model
      ! [GR-SOIL 2026-05-24] Itnumb retired — orphan stub; now state%soilwater%Itnumb.
      ! [GR-SOIL 2026-05-24] fldumpconvcrit retired — orphan stub; now state%cfg%simulation%numerical%dump_convergence_diagnostics.
      ! [GR-SOIL 2026-05-24] flwarn_hc + iwarn_hc retired — orphan stubs; now state%soilwater.
      integer :: dev_cmb  !! Mass balance deviation file unit (previously SAVE in checkmassbal)
      ! [GR-SOIL 2026-05-24] swcaprise retired — orphan stub; now state%cfg%simulation%numerical%swcaprise.
      ! [GR-SOIL 2026-05-24] botcom retired — see state%mesh%botcom
      integer :: dra  !! Internal number of drainage input file *.DRA
      ! [GR-DRA 2026-05-23] dramet retired — see state%drainage%dramet
      integer :: inc  !! Internal number of output file *.INC with incremental water balance data
      ! [GR-DRA 2026-05-23] ipos retired — see state%drainage%ipos
      ! [GR-SOIL 2026-05-24] isoillay retired — read inline from config%soil%isoillay
      ! [GR-SOIL 2026-05-24] ncomp retired — read inline from config%soil%ncomp
      ! [GR-SOIL 2026-05-24] nhead retired — read via state%cfg%soil%initial%z_init size.
      ! [GR-SOIL 2026-05-24] nod1lay retired — see state%mesh%nod1lay
      integer :: nrstaring  !! Number of soil type [1..18] according to Staring series (Wosten et al., 2001)
      ! [GR-SOIL 2026-05-24] nsublay retired — derived inline from config%soil%sublay
      ! [GR-SOIL 2026-05-24] numbit retired — orphan stub; now state%soilwater%numbit.
      ! [GR-SOIL 2026-05-24] numlay retired — see state%mesh%numlay
      integer :: numnodnew  !! Number of desired nodes for soil water quality models
      ! [GR-SOIL 2026-05-24] numtab/numtablay retired — orphan stubs; swsophy=1 dormant.
      integer :: rot  !! Internal number of output file *.ROT with microscopic root water extraction data
      integer :: sw2  !! Switch for prescribed bottom flux: 1 = sine function; 2 = table
      integer :: sw3  !! Switch for prescribed hydraulic head of deep aquifer: 1 = sine function; 2 = table
      ! [GR-SOIL 2026-05-24] sw4 retired — orphan stub; read via state%cfg%bottom_boundary%sw4.
      character(len=1024) :: InList_csv  !! character string with comma-separated list of variables for CSV output
      character(len=1024) :: InList_csv_tz  !! character string with comma-separated list of variables for CSV output
      real(real64), allocatable :: tz_z1_z2(:)  !! Depth range for time-depth CSV output (default: top soil profile, bottom soil profile)
      ! [GR-SOIL 2026-05-24] swbotb3Impl retired — orphan stub; read via state%cfg%bottom_boundary%swbotb3impl.
      ! [GR-BND 2026-05-23] SwBotb3ResVert retired — see state%soilwater%swbotb3resvert
      integer :: swcfbs  !! Switch for use of coefficient CFBS to convert potential ET into potential E: 0 = no; 1 = yes
      integer :: swdiscrvert  !! Switch to convert vertical discretization for soil water quality models: 0 = no; 1 = yes
      ! [GR-DRA 2026-05-23] swdislay/swtopdislay retired — see state%drainage%X
      integer :: swfrost  !! Switch for reduction of hydraulic conductivity in case of frost: 0 = no; 1 = yes
      ! [GR-SOIL 2026-05-24] swhyst retired — orphan stub; read via state%cfg%soil%swhyst.
      ! [GR-SOL 2026-05-24] swinco retired — see state%soilwater%swinco
      ! [GR-DRA 2026-05-23] swliminf retired — see state%drainage%swliminf
      ! [GR-BND 2026-05-23] swpondmx retired — see state%surfacewater%swpondmx
      integer :: swqhbot  !! Switch for flux-groundwater level relationship: 1 = exponential function; 2 = tabular function
      integer :: swcofqhc  !! Switch for additional flux added to exponential flux-groundwater level relationship: 0 = no, 1 = yes
      ! [GR-ATM 2026-05-23] swredu retired — see state%atmosphere%swredu
      ! [GR-SOIL 2026-05-24] ientrytab/ientrytablay retired — orphan stubs; swsophy=1 dormant.
      integer :: swtopsub  !! Switch for topsoil or subsoil: 1 = topsoil, 2 = subsoil
      real(real64) :: aqamp  !! Amplitude of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(real64) :: aqave  !! Average hydraulic head in deep aquifer (L)
      real(real64) :: aqper  !! Period of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(real64) :: aqtmax  !! Time with maximum hydraulic head in deep aquifer (T)
      ! [GR-DRA 2026-05-23] basegw retired — see state%drainage%basegw
      ! [GR-SOL 2026-05-24] bdens retired — see state%soilwater%bdens
      real(real64), allocatable :: c_top(:)  !! Oxygen concentration at top of compartment(kg/m3)
      ! [GR-CROP 2026-05-25] o2_d_soil_term1/term2/gfp100/capac_term/nmin1/mplus1/
      ! ini_stress retired — see state%crop%oxygen (fixed-size macp arrays + flag).
      real(real64) :: cofqha  !! Coefficient A in exponential relationship between drainage flux and groundwater level (L/T)
      real(real64) :: cofqhb  !! Coefficient B in exponential relationship between drainage flux and groundwater level (/T)
      real(real64) :: cofqhc  !! Coefficient C (flux) in exponential relationship between drainage flux and groundwater level (L/T)
      ! [GR-ATM 2026-05-23] cofred retired — see state%atmosphere%cofred
      real(real64) :: CritDevMasBal  !! Maximum error in water balance (L)
      ! [GR-CROP 2026-05-25] CriterHr retired — JvL-only (swdrought=2);
      ! body extracted to src/crop/dormant/jongvanlier.f90.
      real(real64), allocatable :: drares(:)  !! Array with drainage resistance (T) for each drainage level
      real(real64), allocatable :: dznew(:)  !! Desired thickness of compartments for soil water quality models (L)
      ! [GR-DRA 2026-05-23] entres retired — see state%drainage%entres
      ! [GR-DRA 2026-05-23] ftopdislay retired — see state%drainage%ftopdislay
      ! [GR-DRA 2026-05-23] geofac retired — see state%drainage%geofac
      ! [GR-SOIL 2026-05-24] gwli retired — orphan stub; read via state%cfg%soil%gwli.
      ! [GR-BND 2026-05-23] gwltab retired — see state%soilwater%gwltab
      ! [GR-SOIL 2026-05-24] h_enpr retired — orphan stub; lives on state%soilwater%vg_params(:)%h_enpr and config%soil%hydraulics%h_enpr(:).
      ! [GR-BND 2026-05-23] haqtab retired — see state%soilwater%haqtab
      ! [GR-BND 2026-05-23] hbotab retired — see state%soilwater%hbotab
      ! [GR-SOIL 2026-05-24] hcomp retired — derived inline from config%soil%hsublay/ncomp
      real(real64) :: hdrain  !! Mean drainage level (L) to derive regional average groundwater level for bottom boundary condition
      ! [GR-SOIL 2026-05-24] hplate retired — orphan stub; read via state%cfg%bottom_boundary%hplate.
      ! [GR-SOIL 2026-05-24] hsublay retired — read inline from config%soil%hsublay
      real(real64), allocatable :: infres(:)  !! Array with infiltration resistance (T) for each drainage level
      ! [GR-SOL 2026-05-24] inpola/inpolb retired — see state%mesh%{inpola,inpolb}
      real(real64) :: iqinfmax  !! [SS-GR-FINAL D1] qinfmax retired — 0 consumers
      real(real64) :: issnowbeg  !! Amount of snow in soil water equivalent (L) at start of current intermediate period [snow.f90, waterbalance.f90]
      ! [GR-DRA 2026-05-23] khbot retired — see state%drainage%khbot
      ! [GR-DRA 2026-05-23] khtop retired — see state%drainage%khtop
      ! [GR-CROP 2026-05-25] Kroot retired — JvL-only (swdrought=2);
      ! body extracted to src/crop/dormant/jongvanlier.f90.
      ! [GR-SOIL 2026-05-24] ksatthr retired — orphan stub; threshold-Ksat path not ported.
      ! [GR-CROP 2026-05-25] kstem retired — JvL-only (swdrought=2);
      ! body extracted to src/crop/dormant/jongvanlier.f90.
      ! [GR-DRA 2026-05-23] kvbot retired — see state%drainage%kvbot
      ! [GR-DRA 2026-05-23] kvtop retired — see state%drainage%kvtop
      ! [GR-CROP 2026-05-25] OxygenIntercept/OxygenSlope retired —
      ! consumed only by OxygenReproFunction, now dormant
      ! (src/crop/dormant/oxygenrepro.f90); zero live readers/writers.
      ! [GR-CROP 2026-05-25] paramvg retired — tillage now mutates state%soilwater%vg_params_layer(:).
      ! [GR-BND 2026-05-23] pondmxtab retired — see state%surfacewater%pondmxtab
      ! [GR-BND 2026-05-23] qbotab retired — see state%soilwater%qbotab
      ! [GR-SOIL 2026-05-24] qdraincomp retired — orphan stub.
      ! [GR-DRA 2026-05-23] qdrtab retired — see state%drainage%qdrtab
      ! [GR-SOIL 2026-05-24] qimmob retired — orphan stub, fingered-flow flux retired-zero in waterbalance.f90.
      ! [GR-SOIL 2026-05-24] qssdi + qssdisum migrated to state%soilwater (orphan stubs).
      ! [GR-CROP 2026-05-25] dt_SSDI_event migrated to state%crop%irrigation%dt_SSDI_event.
      ! [GR-CROP 2026-05-25] SSDI persistent state migrated to state%crop%irrigation
      ! (16 fields: swssdi_irr, nod_ssdi_irr, ssdi_schedule_irr, ssdi_sched_type_irr,
      !  nod_ssdi_sensor_irr, ssdi_threshold_irr, ssdi_threshold_z_irr, ssdi_amount_irr,
      !  ssdi_appl_rate_irr, sw_interval_irr, days_interval_irr, days_counter_irr,
      !  nirri_ssdi_irr, ssdi_date_irr, ssdi_rate_f_irr, ssdi_amount_f_irr).
      ! [GR-SOIL 2026-05-24] relsatthr retired — orphan stub; threshold-Ksat path not ported.
      ! [GR-SOIL 2026-05-24] rimlay retired — orphan stub; read via state%cfg%bottom_boundary%rimlay.
      ! [GR-ATM 2026-05-23] rsigni retired — see state%atmosphere%rsigni
      ! [GR-ATM 2026-05-23] rsoil retired — see state%atmosphere%rsoil
      ! [GR-ATM 2026-05-23] swuseCN retired — see state%atmosphere%swusecn
      ! [GR-CROP 2026-05-25] Rxylem retired — JvL-only (swdrought=2);
      ! body extracted to src/crop/dormant/jongvanlier.f90.
      ! [GR-BND 2026-05-23] runonarr retired — see state%soilwater%runonarr
      real(real64) :: shape  !! Shape factor: ratio between the mean and the maximum groundwater level elevation above the drainage base (-)
      real(real64) :: sinamp  !! Amplitude of prescribed bottom flux (L/T) in case of sine function
      real(real64) :: sinave  !! Average value of prescribed bottom flux (L/T) in case of sine function
      real(real64) :: sinmax  !! Time of the year with maximum bottom flux in case of prescribed sine function
      ! [GR-SOIL 2026-05-24] sptab/sptablay retired — orphan stubs; swsophy=1 dormant.
      ! [GR-CROP 2026-05-25] StepHr retired — JvL-only (swdrought=2);
      ! body extracted to src/crop/dormant/jongvanlier.f90.
      ! [GR-SOIL 2026-05-24] tau retired — orphan stub; read via state%cfg%soil%tau.
      ! [GR-CROP 2026-05-25] twilt retired — see state%crop%common%twilt
      ! [GR-SOIL 2026-05-24] zi retired — orphan stub; read via state%cfg%soil%initial%z_init.
      ! [GR-DRA 2026-05-23] zintf retired — see state%drainage%zintf
      logical :: FlHydrLift  !! Flag indicating release of water from root to soil is allowed
      ! [GR-BND 2026-05-23] flrunon retired — see state%soilwater%flrunon
      character(len=16) :: drfil  !! Name of drainage input file
      character(len=80) :: pathdrain  !! Path to folder with drainage input files
      integer :: swbotbhea  !! Switch for bottom boundary condition: 1 = heat flux is zero; 2 = prescribed temperature
      integer :: swtopbhea  !! Switch for top boundary condition: 1 = use air temperatures; 2 = read measured surface temperatures
      integer :: swcalt  !! Switch for method of soil water heat flow simulation: 1 = analytical method; 2 = numerical method
      integer :: swhea  !! Switch for simulation of soil heat flow: 0 = no; 1 = yes
      integer :: tem  !! Internal number of output file *.TEM with soil temperatures
      real(real64) :: ddamp  !! Damping depth (L) of temperature wave in soil
      real(real64) :: tampli  !! Amplitude of prescribed annual temperature wave (�C) at soil surface
      real(real64), allocatable :: tembtab(:)  !! Array with specified bottom temperature (�C) as function of time (T)
      real(real64), allocatable :: temtoptab(:)  !! Array with specified soil surface temperature (�C) as function of time (T)
      real(real64) :: tfroststa  !! Soil temperature (�C) where reduction of water fluxes starts
      real(real64) :: tfrostend  !! Soil temperature (�C) where reduction of water fluxes ends
      real(real64) :: timref  !! Time in the year (T) with top of prescribed sine temperature wave
      real(real64) :: tmean  !! Prescribed mean annual temperature (�C) at soil surface
      real(real64), allocatable :: tsoil(:)  !! Config-staging: nheat initial temperature profile entries (�C) for each compartment
      real(real64), allocatable :: zh(:)  !! Array with soil depths (L) used to specify initial soil temperatures
      integer :: snw  !! Internal number of output file *.SNW with snow pack data
      integer :: swsnow  !! Switch for simulation of snow accumulation and melt: 0 = no; 1 = yes
      integer :: swsublim  !! Switch for suppressing simulation of sublimation of snow: 1 = suppress !! Adaptation 3 for PEARL-MACRO
      real(real64) :: snowcoef  !! Snow melt factor (-)
      ! [GR-SOL 2026-05-24] nconc retired — see state%solute%nconc
      ! [GR-SOL 2026-05-24] swbr retired — see state%solute%swbr
      ! [GR-SOL 2026-05-24] swbotbc retired — see state%solute%swbotbc
      integer :: swsolu  !! Switch for simulation of solute transport: 0 = no; 1 = yes
      ! [GR-SOL 2026-05-24] bexp retired — see state%solute%bexp
      ! [GR-SOL 2026-05-24] cdrain retired — see state%solute%cdrain
      ! [GR-SOL 2026-05-24] cirr retired — see state%solute%cirr
      ! [GR-SOL 2026-05-24] cml retired — see state%solute%cml_init / state%solute%cml
      real(real64), allocatable :: cmsy(:)  !! Array with dissolved + adsorbed solute concentration (M/L3 soil volume) in mobile region
      ! [GR-SOL 2026-05-24] cpre retired — see state%solute%cpre
      ! [GR-SOL 2026-05-24] cref retired — see state%solute%cref
      ! [GR-SOL 2026-05-24] cseeptab retired — see state%solute%cseeptab
      ! [GR-SOL 2026-05-24] daquif retired — see state%solute%daquif
      ! [GR-SOL 2026-05-24] ddif retired — see state%solute%ddif
      ! [GR-SOL 2026-05-24] decpot retired — see state%solute%decpot
      ! [GR-SOL 2026-05-24] decsat retired — see state%solute%decsat
      ! [GR-SOL 2026-05-24] dtsolu retired — see state%solute%dtsolu
      ! [GR-SOL 2026-05-24] fdepth retired — see state%solute%fdepth
      ! [GR-SOL 2026-05-24] frexp retired — see state%solute%frexp
      ! [GR-SOL 2026-05-24] gampar retired — see state%solute%gampar
      ! [GR-SOL 2026-05-24] isqbot retired — see state%solute%isqbot
      ! [GR-SOL 2026-05-24] isqtop retired — see state%solute%isqtop
      ! [GR-SOL 2026-05-24] kf retired — see state%solute%kf
      ! [GR-SOL 2026-05-24] kfsat retired — see state%solute%kfsat
      ! [GR-SOL 2026-05-24] ldis retired — see state%solute%ldis
      ! [GR-SOL 2026-05-24] poros retired — see state%solute%poros
      ! [GR-SOL 2026-05-24] rottot retired — see state%solute%rottot
      ! [GR-SOL 2026-05-24] rtheta retired — see state%solute%rtheta
      ! [GR-SOL 2026-05-24] samini retired — see state%solute%samini
      ! [GR-SOL 2026-05-24] sqdra retired — see state%solute%sqdra
      ! [GR-SOL 2026-05-24] tscf retired — see state%solute%tscf
      ! [GR-SOL 2026-05-24] zc retired — see state%solute%zc_init
      real(real64) :: Ageirr  !! Age of irrigation water (d)
      real(real64) :: Agedrain  !! Age of drainage water (d)
      real(real64) :: Agepre  !! Age of precipitation (d)
      real(real64) :: Agepond  !! Age of ponding water (d)
      real(real64) :: Agepondm1  !! Age of ponding water previous timestep (d)
      real(real64) :: icAgetopupw  !! Incremental age leaving top compartment upward (d)
      real(real64) :: icAgetopdwn  !! Incremental age entering top compartment downward (d)
      real(real64) :: Z_Tp  !! [retired-zero] kept: ArMpTp/ArMpSs gating
      real(real64) :: CritUndSatVol  !! [retired-zero] kept: waterbalance watertable() arg
      real(real64) :: ArMpTp  !! [retired-zero] kept: ArMpSs assignment
      ! [GR-SOIL 2026-05-24] cQMpLatSs retired — orphan stub, ADR 0040 macropore.
      real(real64), allocatable :: DiPoCp(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64) :: iQMpOutDrRap  !! [retired-zero] kept: swap_csv/swapoutput DRAINAGE accumulator
      real(real64), allocatable :: IAvFrMpWlWtDm1(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: IAvFrMpWlWtDm2(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: iQExcMtxDm1Cp(:)  !! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(real64), allocatable :: iQExcMtxDm2Cp(:)  !! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(real64), allocatable :: iQOutDrRapCp(:)  !! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(real64), allocatable :: VlMpStDm1(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: VlMpStDm2(:)  !! [retired-zero] kept: soilgrid refinement
      ! [GR-SOIL 2026-05-24] QExcMpMtx + QMaPo retired — orphan stubs, ADR 0040 macropore terms inlined zero in waterbalance.f90.
      ! [GR-DRA 2026-05-23] QRapDra retired — see state%drainage%QRapDra
      ! [GR-DRA 2026-05-23] NumLevRapDra retired — see state%drainage%NumLevRapDra
      logical :: FlDecMpRat  !! [retired-zero] kept: soilhydraulics convergence sentinel
      ! [GR-DRA 2026-05-23] intwl retired — see state%surfacewater%intwl
      ! [GR-DRA 2026-05-23] nowltab retired — see state%drainage%nowltab
      ! [GR-DRA 2026-05-23] impend retired — see state%surfacewater%impend
      ! [GR-DRA 2026-05-23] wlsman retired — see state%surfacewater%wlsman
      ! [GR-DRA 2026-05-23] wlstab retired — see state%surfacewater%wlstab
      real(real64), allocatable :: owltab(:,:)  !! real(8) qdrd                  !! Moved to drainage_state_t%qdrd (ADR 0031)
      ! [GR-CROP 2026-05-25] flCropPrep retired — see state%crop%common%flCropPrep
      ! [GR-CROP 2026-05-25] zPrep/hPrep/MaxPrepDelay/dhPrep retired — see cropwofost_config_t%preparation.
      ! [GR-CROP 2026-05-25] PrepDelay retired — see state%crop%common%PrepDelay
      ! [GR-CROP 2026-05-25] flCropSow retired — see state%crop%common%flCropSow
      ! [GR-CROP 2026-05-25] zSow/hSow/zTempSow/MaxSowDelay/TempSow/dhSow/dtempSow retired — see cropwofost_config_t%sowing.
      ! [GR-CROP 2026-05-25] SowDelay retired — see state%crop%common%SowDelay
      ! [GR-CROP 2026-05-25] flCropGerm retired — see state%crop%common%flCropGerm
      ! [GR-CROP 2026-05-25] zgerm retired — see cropwofost_config_t%germination.
   end type legacy_state_t

end module legacy_state_mod
