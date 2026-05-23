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
      real(real64) :: o2_w_root  !! Dry weight per root length (kg/m)
      real(real64) :: o2_w_root_z0  !! Root weight at depth
      real(real64) :: o2_soil_temp  !! Soil temperature (K)
      real(real64) :: o2_sat_water_cont  !! Saturated water content
      real(real64) :: o2_gas_filled_porosity  !! Gas-filled porosity
      real(real64) :: o2_d_o2inwater  !! O2 diffusion in water
      real(real64) :: o2_d_root  !! Diffusion in root
      real(real64) :: o2_d_soil  !! Soil diffusion
      real(real64) :: o2_perc_org_mat  !! Organic matter percentage
      real(real64) :: o2_soil_density  !! Soil density (kg/m3)
      real(real64) :: o2_depth  !! Compartment thickness (m)
      real(real64) :: o2_shape_factor_microbialr  !! Shape factor microbial resp
      real(real64) :: o2_root_radius  !! Root radius (m)
      real(real64) :: o2_r_microbial_z0  !! Microbial respiration rate
      real(real64) :: o2_waterfilm_thickness  !! Water film thickness
      real(real64) :: o2_bunsencoeff  !! Bunsen coefficient
      real(real64) :: o2_c_min_micro  !! Min O2 for microbial resp
      real(real64) :: o2_c_macro  !! Macropore O2 conc
      real(real64) :: o2_ctopnode  !! Top node O2 conc
      real(real64), allocatable :: outdat(:)  !! Array with output dates for water and solute balances
      real(real64), allocatable :: outdatint(:)  !! Array with intermediate output dates
      logical :: flSwapShared  !! Flag to indicate the shared simultaneous simulation with other applications
      character(len=16) :: outfil  !! Name of output file
      character(len=80) :: pathwork  !! Path to work directory
      character(len=80) :: project  !! Name of project
      integer :: daynrfirst  !! First calendar day number for which meteorological data is available in current year
      integer :: daynrlast  !! Last calendar day number for which meteorological data is available in current year
      ! [GR-ATM 2026-05-23] detrecord/irectotal retired — see state%atmosphere%{detrecord,irectotal}
      integer :: nofd  !! number of days for running average Tmin (-)
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
      integer :: irrigevent  !! Switch: 0 = no irrigation; 1 = fixed irrigation event; 2 = scheduled irrigation event
      integer, allocatable :: irtype(:)  !! Type of fixed irrigation: 0 = sprinkling irrigation; 1 = surface irrigation
      integer :: isua  !! Switch for type of irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer :: isuas  !! Switch for type of scheduled irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer :: nirri  !! Number of irrigation event
      integer :: swirfix  !! Switch for fixed irrigation: 0 = no applications prescribed; 1 = applications are prescribed
      integer :: swcirrthres  !! Switch to allow over irrigation when a conc-threshold is exceeded: 0 = no; 1 = yes/allowed
      real(real64) :: cirrs  !! Solute concentration of irrigation water (M/L3)
      real(real64) :: cirrthres  !! Threshold value (M/L3) indicating the concentration that initiates over irrigation
      real(real64) :: dcrit  !! Depth (L) of sensor for soil water pressure head or water content
      real(real64), allocatable :: ditab(:)  !! Array with amount of under- or over-irrigation (L) as function of crop development stage
      real(real64), allocatable :: dwatab(:)  !! Array with maximum amounts of water depleted as function of crop development stage
      real(real64), allocatable :: fidtab(:)  !! Array with prescribed fixed irrigation depth (L) as function of crop development stage
      real(real64), allocatable :: hcritab(:)  !! Array with minimum soil water pressure heads (L) as function of crop development stage
      real(real64), allocatable :: irconc(:)  !! Array with irrigation concentrations (M/L3) in case of fixed irrigation
      real(real64), allocatable :: irdate(:)  !! Array with fixed irrigation dates
      real(real64), allocatable :: irdepth(:)  !! Array with fixed irrigation depths (L)
      ! [GR-ATM 2026-05-23] nird retired — see state%atmosphere%nird
      real(real64) :: perirrsurp  !! percentage (-) of the scheduled irrigation depths that may be over irrigated
      real(real64) :: raithreshold  !! Threshold value (L) indicating the amount of rainfall which is substracted from scheduled irrigation depths
      real(real64), allocatable :: rawtab(:)  !! Array with minimum of readily available water as function of crop development stage
      real(real64), allocatable :: tawtab(:)  !! Array with minimum of totally available water as function of crop development stage
      real(real64), allocatable :: tcritab(:)  !! Array with minimum volumetric soil water contents as function of crop development stage
      real(real64) :: tstairrig  !! Date after which scheduled irrigation is allowed
      real(real64) :: tendirrig  !! Date after which scheduled irrigation is NOT allowed
      real(real64), allocatable :: treltab(:)  !! Array with minimum of ratio actual/potential transpiration as function of crop development stage
      integer :: dayfix  !! days since last irrigation event
      integer :: tsumtime  !! time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      real(real64) :: tsumtemp  !! temperature limit to initiate grass growth  [0.0..20.0 grC, R]
      real(real64) :: tsumdepth  !! depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
      logical :: flCropCalendar  !! Flag indicating that crop season is active (but currently might be bare or cropped)
      logical :: flCropEmergence  !! Flag indicating period from crop emergence until harvest
      logical :: flCropHarvest  !! Flag indicating period from crop harvest until the end of crop season
      logical :: flCropReadFile  !! Flag indicating reading of input.crp
      logical :: flCropOpenFile  !! Flag indicating to create output.crp
      logical :: flCropOutput  !! Flag indicating writing of output.crp
      ! [GR-ATM 2026-05-23] croptype retired — see state%crop%common%croptype
      integer :: swcrp  !! Switch for output file *.CRP with daily crop output: 0 = no; 1 = yes
      integer :: crp  !! Internal number of crop output file *.CRP
      integer :: daycrop  !! Number of days that a crop exists
      integer :: icrop  !! Current crop number
      integer :: idsl  !! Switch for crop development before anthesis: 0 = depends on temperature;
      integer :: noddrz  !! Compartment number at bottom root zone (-)
      ! [GR-ATM 2026-05-23] atmtr retired — see state%atmosphere%atmtr
      real(real64) :: agerm  !! Coefficient a  of germination
      real(real64) :: cgerm  !! Coefficient c  of germination
      real(real64) :: bgerm  !! Coefficient b  of germination
      integer :: swpotrelmf  !! Calculation of potential yield
      real(real64), allocatable :: avevaptb(:)  !! Gash interception model: average evaporation intensity during shower (-) as function of time (T)
      real(real64), allocatable :: avprectb(:)  !! Gash interception model: average rainfall intensity (-) as function of time (T)
      real(real64) :: c_mroot  !! Maintenance coefficient of root [0.0..1.0 kg O2/kg/d, R]
      real(real64), allocatable :: cropend(:)  !! Array with crop end dates
      real(real64), allocatable :: cropstart(:)  !! Array with crop start dates
      ! [GR-ATM 2026-05-23] difpp retired — see state%atmosphere%difpp
      real(real64) :: dlc  !! Shortest day length (T) for any crop development
      real(real64) :: dlo  !! Minimum day length (T) for optimal crop development
      ! [GR-ATM 2026-05-23] dsinbe retired — see state%atmosphere%dsinbe
      real(real64) :: f_senes  !! Reduction factor for senescence, used for maintenance respiration [0..1.0 -, R]
      real(real64) :: gasst  !! Total gross assimilation for actual crop (kg/ha)
      real(real64) :: gasstpot  !! Total gross assimilation for potential crop (kg/ha)
      ! [GR-ATM 2026-05-23] gc retired — see state%crop%common%gc
      real(real64) :: hdrygerm  !! Criterium Hdry of germination
      real(real64) :: hwetgerm  !! Criterium Hwet of germination
      real(real64) :: max_resp_factor  !! Ratio root total respiration / maintenance respiration [1..5.0 -, R]
      real(real64) :: mrest  !! Total maintenance respiration for actual crop (kg/ha)
      real(real64) :: mrestpot  !! Total maintenance respiration for potential crop (kg/ha)
      real(real64), allocatable :: mrftb(:)  !! Array with ratio root total respiration / maintenance respiration as function of DVS (kg/m3)
      real(real64), allocatable :: pfreetb(:)  !! Gash interception model: free throughfall coefficient (-) as function of time (T)
      real(real64), allocatable :: pstemtb(:)  !! Gash interception model: stem flow coefficient (-) as function of time (T)
      ! [GR-ATM 2026-05-23] siccaptb retired — see state%crop%common%siccaptb
      real(real64) :: siccaplai  !! interception storage per unit of LAI (cm/LAI)
      real(real64) :: q10_microbial  !! Relative increase in microbial respiration at temperature increase of 10 �C [1.0..4.0 -, R]
      real(real64) :: q10_root  !! Relative increase in root respiration at temperature increase of 10 �C [1.0..4.0 -, R]
      real(real64) :: reltr  !! relative transpiration factor that reduces crop growth (-)
      real(real64) :: rid  !! Real day number of detailed grass crop (d)
      real(real64) :: rootcoefa  !! Defines relative distance at which mean soil water content occurs between roots
      real(real64) :: rooteff  !! Root system efficiency factor [0..1.0 -, R]
      real(real64) :: rootradius  !! Root radius drought stress (cm)
      real(real64), allocatable :: scanopytb(:)  !! Gash interception model: storage capacity of canopy (-) as function of time (T)
      real(real64) :: shape_factor_rootr  !! Shape factor for exponential decrease of root respiration rate with depth [0..1.0 -, R]
      real(real64) :: specific_resp_humus  !! Respiration rate of humus at 25 �C [0.0..1.0 kg O2/kg C/d, R]
      real(real64) :: tadw  !! Dry weight of plant minus roots of actual growth (kg/ha)
      real(real64) :: tadwpot  !! Dry weight of plant minus roots of potential growth (kg/ha)
      real(real64) :: tsumemeopt  !! Temperature sum for crop emergence under optimal conditions
      real(real64) :: tsumgerm  !! Temperature sum during germination
      real(real64) :: TBASEM  !! Lower threshold temp. for emergence (C)
      real(real64) :: TEFFMX  !! max. eff. temp. for emergence (C)
      real(real64) :: w_root_ss  !! Dry weight of roots at soil surface [0.0..10.0 kg/m3, R]
      real(real64) :: wiltpoint  !! Minimum pressure head at interface soil-root (cm)
      real(real64), allocatable :: wrtb(:)  !! Array with dry weight of root at soil surface as function of DVS (kg/m3)
      real(real64) :: wrtmin  !! Minimum dry weight of plant root at relative depth (1% of the initial value)
      real(real64) :: gwrt  !! Growth of dry weight of plant root (kg/ha)
      logical :: flanthesis  !! Flag indicating anthesis stage of a crop
      logical :: flHarvestDay  !! Flag indicating that current day is harvest day
      character(len=40), allocatable :: cropfil(:)  !! Array with names of crop files
      character(len=80) :: pathcrop  !! Path to folder with crop input files
      real(real64) :: rdmax  !! Maximum rooting depth in soil profile (L)
      ! [GR-ATM 2026-05-23] flCO2 retired — see state%atmosphere%flco2
      real(real64), allocatable :: co2amaxtb(:)  !! table with factors to correct AMAX for CO2
      real(real64), allocatable :: co2efftb(:)  !! table with factors to correct EFF for CO2
      real(real64), allocatable :: co2tratb(:)  !! table with factors to correct TRA for CO2
      integer, allocatable :: co2year(:)  !! table with years for which CO2 concentrations are given
      real(real64), allocatable :: co2ppm(:)  !! table with CO2 concentrations (ppm), for each year in co2year
      real(real64) :: verndvs  !! critical development stage after which the effect of vernalisation is halted [-]
      real(real64) :: vernsat  !! saturated vernalisation requirement [d]
      real(real64) :: vernbase  !! base vernalisation requirement [d]
      real(real64), allocatable :: vernrtb(:)  !! table with rate of vernalisation as function of tav [days/degrees]
      integer :: swbulb  !! switch to enable simulation of bulb crops (-)
      real(real64) :: drbl  !! Death rate of actual bulb (kg/ha)
      real(real64) :: drblpot  !! Death rate of potential bulb (kg/ha)
      real(real64) :: fbl  !! Dry weight fraction partitioned to flowers (-)
      real(real64) :: pld
      real(real64) :: remoc  !! [GR-CROP-DVS] plwt retired — see state%crop%wofost%plwt
      logical :: flCropNut  !! Flag indicating simulation of crop nutrient stress
      real(real64), allocatable :: nmxlv(:)
      integer :: ilnmxl
      real(real64) :: fstr  !! [SS-GR-FINAL D1] amFERT retired — 0 consumers
      integer :: till_swtill  !! Switch: 0=no tillage, 1=tillage
      integer :: till_i_n_model  !! Switch for n-parameter treatment (1-3)
      integer :: till_iRedist  !! Redistribution type after MvG change
      integer :: till_Ntill  !! Number of tabulated tillage events
      integer :: till_Ntypes  !! Number of tillage types
      real(real64) :: till_Max_Z_tillage  !! Max possible depth of tillage (cm)
      integer, allocatable :: iHWCKmodel(:)  !! indicator what type of water retention and hydraulic conductivity model is used (per soil layer)
      logical, allocatable :: BiModal(:)  !! logical indicating whether chosen model is bi-modal or not
      logical, allocatable :: NoVap(:)  !! logical indicating that NO vapour flow is to be considered in PDI K-model
      integer, allocatable :: Itnumb(:,:)  !! Iteration number statistics [soilhydraulics.f90, timecontrol_mod.f90]
      logical :: fldumpconvcrit  !! flag to generate additional output about convergence-warnings from subr Headcalc
      logical :: flwarn_hc  !! Headcalc warning flag (previously SAVE variable)
      integer :: iwarn_hc  !! Headcalc warning counter (previously SAVE variable)
      integer :: dev_cmb  !! Mass balance deviation file unit (previously SAVE in checkmassbal)
      logical :: swcaprise  !! flag to minimize cap.rise to rootzone (for experts only)
      integer, allocatable :: botcom(:)  !! Array with number of bottom compartments in each soil layer
      integer :: dra  !! Internal number of drainage input file *.DRA
      integer :: dramet  !! Switch for lateral drainage: 1 = table of flux - groundwater level; 2 = Hooghoudt or Ernst;
      integer :: inc  !! Internal number of output file *.INC with incremental water balance data
      integer :: ipos  !! Switch for position of drain (see *.DRA input file for overview)
      integer, allocatable :: isoillay(:)  !! Number of soil layer, starting with 1 at the soil surface
      integer, allocatable :: ncomp(:)  !! Array with number of compartments in each sublayer
      integer :: nhead  !! Number of initial soil water pressure heads as provided in the input
      integer, allocatable :: nod1lay(:)  !! node nr of first node of each soil layer (from top to bottom)
      integer :: nrstaring  !! Number of soil type [1..18] according to Staring series (Wosten et al., 2001)
      integer :: nsublay  !! Number of sublayers in the soil profile
      integer :: numbit  !! Iteration number for solving Richards equation
      integer :: numlay  !! Number of (physical) soil layers
      integer :: numnodnew  !! Number of desired nodes for soil water quality models
      integer, allocatable :: numtab(:)  !! Number of table entries of soil physical values for each model compartment
      integer, allocatable :: numtablay(:)  !! Number of table entries of soil physical values for each soil layer
      integer :: rot  !! Internal number of output file *.ROT with microscopic root water extraction data
      integer :: sw2  !! Switch for prescribed bottom flux: 1 = sine function; 2 = table
      integer :: sw3  !! Switch for prescribed hydraulic head of deep aquifer: 1 = sine function; 2 = table
      integer :: sw4  !! Switch for extra groundwater flux as function of time: 0 = no extra flux; 1 = include extra flux
      character(len=1024) :: InList_csv  !! character string with comma-separated list of variables for CSV output
      character(len=1024) :: InList_csv_tz  !! character string with comma-separated list of variables for CSV output
      real(real64), allocatable :: tz_z1_z2(:)  !! Depth range for time-depth CSV output (default: top soil profile, bottom soil profile)
      integer :: swbotb3Impl  !! Switch for implicit solution with lower boundary option 3 (Cauchy): 0 = explicit, 1 = implicit
      ! [GR-BND 2026-05-23] SwBotb3ResVert retired — see state%soilwater%swbotb3resvert
      integer :: swcfbs  !! Switch for use of coefficient CFBS to convert potential ET into potential E: 0 = no; 1 = yes
      integer :: swdiscrvert  !! Switch to convert vertical discretization for soil water quality models: 0 = no; 1 = yes
      integer :: swdislay  !! Switch to distribute drainage flux vertically with a given position of the top of the model discharge layers: 0 = no; 1 = yes
      integer, allocatable :: swtopdislay(:)  !! Switch, for each drainage level, to distribute drainage flux vertically with a given position of the top of the model discharge layers: 0 = no; 1 = yes
      integer :: swfrost  !! Switch for reduction of hydraulic conductivity in case of frost: 0 = no; 1 = yes
      integer :: swhyst  !! Switch for hysteresis of soil moisture retention function: 0 = no; 1 = yes
      integer :: swinco  !! Switch for initial soil moisture condition: 1 = pressure heads; 2 = hydrostatic equilibrium;
      integer :: swliminf  !! Switch for limit of infiltration head to the waterdepth in the channel: 0 = nolimit, 1 = limitation
      ! [GR-BND 2026-05-23] swpondmx retired — see state%surfacewater%swpondmx
      integer :: swqhbot  !! Switch for flux-groundwater level relationship: 1 = exponential function; 2 = tabular function
      integer :: swcofqhc  !! Switch for additional flux added to exponential flux-groundwater level relationship: 0 = no, 1 = yes
      ! [GR-ATM 2026-05-23] swredu retired — see state%atmosphere%swredu
      integer, allocatable :: ientrytab(:,:)  !! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      integer, allocatable :: ientrytablay(:,:)  !! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer
      integer :: swtopsub  !! Switch for topsoil or subsoil: 1 = topsoil, 2 = subsoil
      real(real64) :: aqamp  !! Amplitude of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(real64) :: aqave  !! Average hydraulic head in deep aquifer (L)
      real(real64) :: aqper  !! Period of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(real64) :: aqtmax  !! Time with maximum hydraulic head in deep aquifer (T)
      real(real64) :: basegw  !! Depth of impervious layer (L) for drainage according to Hooghoudt or Ernst
      real(real64), allocatable :: bdens(:)  !! Array with dry bulk density for each soil layer (M/L3)
      real(real64), allocatable :: c_top(:)  !! Oxygen concentration at top of compartment(kg/m3)
      real(real64), allocatable :: o2_d_soil_term1(:)  !! Pre-calculated soil diffusion term1 per node
      real(real64), allocatable :: o2_d_soil_term2(:)  !! Pre-calculated soil diffusion term2 per node
      real(real64), allocatable :: o2_gfp100(:)  !! Gas-filled porosity * 100 per node
      real(real64), allocatable :: o2_capac_term(:)  !! Water capacity term per node
      real(real64), allocatable :: o2_nmin1(:)  !! N-1 per node for VG equation
      real(real64), allocatable :: o2_mplus1(:)  !! M+1 per node for VG equation
      logical :: o2_ini_stress  !! O2 stress initialization flag (initialized to .true. via data statement)
      real(real64) :: cofqha  !! Coefficient A in exponential relationship between drainage flux and groundwater level (L/T)
      real(real64) :: cofqhb  !! Coefficient B in exponential relationship between drainage flux and groundwater level (/T)
      real(real64) :: cofqhc  !! Coefficient C (flux) in exponential relationship between drainage flux and groundwater level (L/T)
      ! [GR-ATM 2026-05-23] cofred retired — see state%atmosphere%cofred
      real(real64) :: CritDevMasBal  !! Maximum error in water balance (L)
      real(real64) :: CriterHr  !! Maximum difference of Hroot between iterations; convergence criterium  (L)
      real(real64), allocatable :: drares(:)  !! Array with drainage resistance (T) for each drainage level
      real(real64), allocatable :: dznew(:)  !! Desired thickness of compartments for soil water quality models (L)
      real(real64) :: entres  !! Drain entry resistance (T)
      real(real64), allocatable :: ftopdislay(:)  !! Array with factor for function to determine depth of top of model discharge layer for each drain level, see also swtopdislay (L)
      real(real64) :: geofac  !! Geometry factor (-) for analytical drainage formula of Ernst
      real(real64) :: gwli  !! Groundwater level (L) at start of simulation
      ! [GR-BND 2026-05-23] gwltab retired — see state%soilwater%gwltab
      real(real64), allocatable :: h_enpr(:)  !! Soil water Entry Pressure head for Modified MualemVanGenuchten curve (L)
      ! [GR-BND 2026-05-23] haqtab retired — see state%soilwater%haqtab
      ! [GR-BND 2026-05-23] hbotab retired — see state%soilwater%hbotab
      real(real64), allocatable :: hcomp(:)  !! Array with prescribed height of numerical compartments (L) for each sublayer
      real(real64) :: hdrain  !! Mean drainage level (L) to derive regional average groundwater level for bottom boundary condition
      real(real64) :: hplate  !! Pressure head of ceramic plate below lysimeter
      real(real64), allocatable :: hsublay(:)  !! Array with prescribed height of sublayers (L)
      real(real64), allocatable :: infres(:)  !! Array with infiltration resistance (T) for each drainage level
      real(real64), allocatable :: inpola(:)  !! Weight for interpolation between current node and upper node
      real(real64), allocatable :: inpolb(:)  !! Weight for interpolation between current node and lower node
      real(real64) :: iqinfmax  !! [SS-GR-FINAL D1] qinfmax retired — 0 consumers
      real(real64) :: issnowbeg  !! Amount of snow in soil water equivalent (L) at start of current intermediate period [snow.f90, waterbalance.f90]
      real(real64) :: khbot  !! Horizontal hydraulic conductivity of bottom layer (L/T)
      real(real64) :: khtop  !! Horizontal hydraulic conductivity of top layer (L/T)
      real(real64) :: Kroot  !! Hydraulic, radial conductivity of root tissue (L/T)
      real(real64), allocatable :: ksatthr(:)  !! Array with saturated hydraulic conductivity (L/T) for each soil layer: to interpolate VG and Ksatexm
      real(real64) :: kstem  !! Conductance in the path from leaf to root xylem (/d)
      real(real64) :: kvbot  !! Vertical hydraulic conductivity of bottom layer (L/T)
      real(real64) :: kvtop  !! Vertical hydraulic conductivity of top layer (L/T)
      real(real64), allocatable :: OxygenIntercept(:)  !! Parameters of reproduction function for oxygen stress according to Bartholomeus
      real(real64), allocatable :: OxygenSlope(:)  !! Parameters of reproduction function for oxygen stress according to Bartholomeus
      real(real64), allocatable :: paramvg(:,:)  !! Array with input values of soil hydraulic parameters according to Mualem - van Genuchten for each soil layer
      ! [GR-BND 2026-05-23] pondmxtab retired — see state%surfacewater%pondmxtab
      real(real64), allocatable :: qbotab(:)  !! Array with specified bottom flux (L/T) as function of time (T)
      real(real64), allocatable :: qdraincomp(:)  !! Total lateral drainage flux (L/T) for each compartment
      real(real64), allocatable :: qdrtab(:)  !! Array with lateral drainage flux (L/T) as function of groundwater level (L)
      real(real64), allocatable :: qimmob(:)  !! Soil water flux between mobile and immobile fraction in case of fingered flow (L/T)
      real(real64) :: qssdisum  !! Total subsurface irrigation flux (L/T)
      real(real64), allocatable :: qssdi(:)  !! Array with water input via subsurface drip irrigation for each compartment (L/T)
      real(real64) :: dt_SSDI_event  !! Length of SSDI irrigation event (T)
      integer :: swssdi_irr  !! Switch: SSDI active (0=no, 1=yes)
      integer, allocatable :: nod_ssdi_irr(:)  !! Upper and lower nodes for SSDI
      integer :: ssdi_schedule_irr  !! Schedule type (0=fixed dates, 1=internal)
      integer :: ssdi_sched_type_irr  !! Internal schedule type (1=Tact/Tpot, 2=h, 3=theta)
      integer :: nod_ssdi_sensor_irr  !! Sensor node (if ssdi_sched_type > 1)
      real(real64) :: ssdi_threshold_irr  !! Threshold value for scheduling
      real(real64) :: ssdi_threshold_z_irr  !! Depth for threshold value
      real(real64) :: ssdi_amount_irr  !! Amount of scheduled irrigation (cm)
      real(real64) :: ssdi_appl_rate_irr  !! Application rate (cm/d)
      integer :: sw_interval_irr  !! Switch for minimum interval
      integer :: days_interval_irr  !! Minimum days between applications
      integer :: days_counter_irr  !! Days since previous application
      integer :: nirri_ssdi_irr  !! SSDI counter/entry point
      real(real64), allocatable :: ssdi_date_irr(:)  !! Fixed irrigation dates
      real(real64), allocatable :: ssdi_rate_f_irr(:)  !! Fixed irrigation rates (cm/d)
      real(real64), allocatable :: ssdi_amount_f_irr(:)  !! Fixed irrigation amounts (cm)
      real(real64), allocatable :: relsatthr(:)  !! Array with relative saturation (-) for each soil layer: to interpolate VG and Ksatexm
      real(real64) :: rimlay  !! Vertical resistance of aquitard (T)
      ! [GR-ATM 2026-05-23] rsigni retired — see state%atmosphere%rsigni
      ! [GR-ATM 2026-05-23] rsoil retired — see state%atmosphere%rsoil
      ! [GR-ATM 2026-05-23] swuseCN retired — see state%atmosphere%swusecn
      real(real64) :: Rxylem  !! Mean radius of xylem tube inside roots (L)
      ! [GR-BND 2026-05-23] runonarr retired — see state%soilwater%runonarr
      real(real64) :: shape  !! Shape factor: ratio between the mean and the maximum groundwater level elevation above the drainage base (-)
      real(real64) :: sinamp  !! Amplitude of prescribed bottom flux (L/T) in case of sine function
      real(real64) :: sinave  !! Average value of prescribed bottom flux (L/T) in case of sine function
      real(real64) :: sinmax  !! Time of the year with maximum bottom flux in case of prescribed sine function
      real(real64), allocatable :: sptab(:,:,:)  !! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      real(real64), allocatable :: sptablay(:,:,:)  !! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer
      real(real64) :: StepHr  !! Maximum difference of Hroot and Hxylem between iterations; convergence criterium  (L)
      real(real64) :: tau  !! Minimum pressure head difference (L) to change from wetting to drying in case of hysteresis
      real(real64), allocatable :: twilt(:)  !! Pressure head of a compartment at wilting point (L)
      real(real64), allocatable :: zi(:)  !! Array with soil depths (L) used to specify initial soil water pressure heads
      real(real64) :: zintf  !! Depth (L) at which fine top layer ends and coarse sub layer starts
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
      integer :: nconc  !! Number of initial solute concentrations as provided in the input
      integer :: swbr  !! Switch to consider mixed reservoir for solute breakthrough in the saturated zone: 0 = no; 1 = yes
      integer :: swbotbc  !! Switch for bottom boundary condition of solute-concentration (see *.SWP input file for overview)
      integer :: swsolu  !! Switch for simulation of solute transport: 0 = no; 1 = yes
      real(real64) :: bexp  !! Exponent in decomposition reduction factor due to dryness (-)
      real(real64) :: cdrain  !! Mean solute concentration in aquifer or drainage system (M/L3 water)  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64) :: cirr  !! Solute concentration (M/L3) in irrigation water
      real(real64), allocatable :: cml(:)  !! Array with solute concentration (M/L3 water) in mobile region
      real(real64), allocatable :: cmsy(:)  !! Array with dissolved + adsorbed solute concentration (M/L3 soil volume) in mobile region
      real(real64) :: cpre  !! Solute concentration (M/L3) in precipitation
      real(real64) :: cref  !! Reference solute concentration (M/L3) for Freundlich adsorption
      real(real64), allocatable :: cseeptab(:)  !! Array with Mean solute concentration in upward seepage water at bottom of profile (M/L3 water) as function of time (T)
      real(real64) :: daquif  !! Thickness of saturated aquifer (L) to calculate solute breakthrough to surface water
      real(real64) :: ddif  !! Molecular diffusion coefficient (L2/T)
      real(real64), allocatable :: decpot(:)  !! Array with Potential decomposition rate (/T) for each soil layer
      real(real64) :: decsat  !! Decomposition rate in aquifer (/T)
      real(real64) :: dtsolu  !! Maximum time step (T) for accurate numerical solution of solute transport equation  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64), allocatable :: fdepth(:)  !! Array with reduction factor for decomposition (-) for each soil layer
      real(real64) :: frexp  !! Array with Freundlich exponent (-) for solute adsorption
      real(real64) :: gampar  !! Reduction factor for decomposition due to low temperatures (/C)
      real(real64) :: isqbot  !! Solute flux at the bottom of the soil column (M/L2/T)  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64) :: isqtop  !! Solute flux through the soil top surface (M/L2/T)  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64), allocatable :: kf(:)  !! Array with Freundlich coefficient (L3/M) for solute adsorption for each soil layer
      real(real64) :: kfsat  !! Linear adsorption coefficient in aquifer (L3/M)
      real(real64), allocatable :: ldis(:)  !! Array with Solute dispersion length (L) for each soil layer
      real(real64) :: poros  !! Porosity of aquifer (-) to calculate solute breakthrough
      real(real64) :: rottot  !! Cumulative amount of solutes (M/L2) extracted by plant roots  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64) :: rtheta  !! Minimum volumetric water content (-) for potential decomposition
      real(real64) :: samini  !! Total amount of solutes (M/L2) in soil profile at start of current balance period  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64) :: sqdra  !! Total amount of solutes (M/L2) transported to drainage canals  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(real64) :: tscf  !! Relative uptake of solutes by roots (-)
      real(real64), allocatable :: zc(:)  !! Array with soil depths (L) used to specify initial solute concentrations
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
      real(real64) :: cQMpLatSs  !! [retired-zero] kept: soilhydraulics zero-write
      real(real64), allocatable :: DiPoCp(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64) :: iQMpOutDrRap  !! [retired-zero] kept: swap_csv/swapoutput DRAINAGE accumulator
      real(real64), allocatable :: IAvFrMpWlWtDm1(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: IAvFrMpWlWtDm2(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: iQExcMtxDm1Cp(:)  !! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(real64), allocatable :: iQExcMtxDm2Cp(:)  !! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(real64), allocatable :: iQOutDrRapCp(:)  !! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(real64), allocatable :: VlMpStDm1(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: VlMpStDm2(:)  !! [retired-zero] kept: soilgrid refinement
      real(real64), allocatable :: QExcMpMtx(:)  !! [retired-zero] kept: waterbalance use clause
      real(real64) :: QMaPo  !! [retired-zero] kept: waterbalance qbot term
      real(real64) :: QRapDra  !! [retired-zero] kept: surfacewater drainage terms
      integer :: NumLevRapDra  !! Number of drainage levels for rapid drainage (drainage feature)
      logical :: FlDecMpRat  !! [retired-zero] kept: soilhydraulics convergence sentinel
      integer, allocatable :: intwl(:)  !! imper removed (surfacewater_state_t%imper)
      integer, allocatable :: nowltab(:)
      real(real64), allocatable :: impend(:)
      real(real64), allocatable :: wlsman(:,:)  !! SS-SWST Phase 2 Task 11 C2: sttab removed — state%surfacewater%sttab owns it.
      real(real64), allocatable :: wlstab(:)  !! sttab(22,2) removed (surfacewater_state_t%sttab)
      real(real64), allocatable :: owltab(:,:)  !! real(8) qdrd                  !! Moved to drainage_state_t%qdrd (ADR 0031)
      logical :: flCropPrep  !! Flag indicating if ploughing opportunity has been realized
      real(real64) :: zPrep  !! z-level for monitoring work-ability for the crop
      real(real64) :: hPrep  !! maximum pressure head during preparation
      integer :: MaxPrepDelay  !! maximum delay of preparation (starting from begin of growing season)
      integer :: PrepDelay  !! delay of preparation
      real(real64) :: dhPrep  !! overshoot of pressure head for work-ability during preparation
      logical :: flCropSow  !! Flag indicating if sowing opportunity has been realized
      real(real64) :: zSow  !! z-level for monitoring work-ability for the crop
      real(real64) :: hSow  !! maximum pressure head during sowing
      real(real64) :: zTempSow  !! z-level for monitoring temperature for sowing
      integer :: MaxSowDelay  !! maximum delay of sowing (starting from begin of growing season)
      integer :: SowDelay  !! delay of delay
      real(real64) :: TempSow  !! temperature for sowing
      real(real64) :: dhSow  !! overshoot of pressure head for work-ability during sowing
      real(real64) :: dtempSow  !! undershoot of temperature for sowing at end of available period
      logical :: flCropGerm  !! Flag indicating if germination has been realized
      real(real64) :: zgerm  !! z-level for monitoring temperature for germination
   end type legacy_state_t

end module legacy_state_mod
