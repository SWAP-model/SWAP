! File VersionID:
!   $Id: variables.f90 378 2018-05-08 13:50:52Z heine003 $
! ----------------------------------------------------------------------
! File Name:    variables.for
! Content:      common variables of modular SWAP code
! Sections:     time & control, meteo, irrigation, crop, soilwater, macropore, surfacewater, heat, snow, solute
      module  variables
      implicit none
      save
      include 'arrays.fi'

! --- time & control variables
      ! ========================================================================
      ! [SS-SWC] retired 2026-05-12 — soil-water core migrated to state%soilwater
      !   See: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md
      !        ADR 0038 (next task S-2.13)
      ! ========================================================================
      ! real(8)   iqredwet_day       ! moved to state%soilwater%iqredwet_day
      ! real(8)   iqreddry_day       ! moved to state%soilwater%iqreddry_day
      ! real(8)   iqredsol_day       ! moved to state%soilwater%iqredsol_day
      ! real(8)   iqredfrs_day       ! moved to state%soilwater%iqredfrs_day
      ! real(8)   iptra_day          ! moved to state%soilwater%iptra_day

      ! ========================================================================
      ! [SS-TC] retired 2026-05-12 — TimeControl runtime fields moved to
      !   state%timecontrol. See:
      !     docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
      !     ADR 0041 (next task TC-15)
      ! ========================================================================
      ! TimeControl persistent state (formerly local SAVE variables in timecontrol.f90)
      ! integer   tc_datea(6)        ! moved to state%timecontrol%datea
      ! integer   tc_nextyear        ! moved to state%timecontrol%nextyear
      ! integer   tc_flprevious      ! moved to state%timecontrol%flprevious
      ! logical   tc_flTnext         ! moved to state%timecontrol%flTnext
      ! real(4)   tc_fsec            ! moved to state%timecontrol%fsec
      ! real(8)   tc_tchange         ! moved to state%timecontrol%tchange
      ! real(8)   tc_dtEvent         ! moved to state%timecontrol%dtEvent
      ! real(8)   tc_tEvent          ! moved to state%timecontrol%tEvent
      ! real(8)   tc_tcumold         ! moved to state%timecontrol%tcumold
      ! real(8)   tc_dtprevious      ! moved to state%timecontrol%dtprevious
      ! real(4)   tc_tmptimestart    ! moved to state%timecontrol%tmptimestart
      ! real(4)   tc_tmptimeend      ! moved to state%timecontrol%tmptimeend
      ! logical   FlOpenFileDev      ! moved to state%timecontrol%flOpenFileDev
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%* (ADR 0041)
      ! integer   nprintday          ! moved to state%timecontrol%nprintday
      ! logical   flprintdt          ! moved to state%timecontrol%flprintdt
      ! logical   flprintshort       ! moved to state%timecontrol%flprintshort
      ! logical   floutputshort      ! moved to state%timecontrol%floutputshort
      ! integer   nprintcount        ! moved to state%timecontrol%nprintcount
      ! integer   cntper             ! moved to state%timecontrol%cntper
      ! integer   daycum             ! moved to state%timecontrol%daycum
      ! integer   daynr              ! moved to state%timecontrol%daynr
      ! integer   imonth             ! moved to state%timecontrol%imonth
      ! integer   ioutdat            ! moved to state%timecontrol%ioutdat
      ! integer   ioutdatint         ! moved to state%timecontrol%ioutdatint
      ! integer   isteps             ! moved to state%timecontrol%isteps
      ! integer   iyear              ! moved to state%timecontrol%iyear
      ! integer   iyearm1            ! moved to state%timecontrol%iyearm1
      integer   logf               ! Internal number of logbook output file *.LOG

      ! Exchange/DLL persistent state (formerly local SAVE in swap.f90 handle_exchange)
      real(8)   ex_tlast           ! [SS-DRV Task 6] dead — handle_exchange retired; only initialize.f90 zeros it. Retirement candidate.
      
      ! Oxygen stress persistent state (moved from O2_pars module and OxygenStress subroutine)
      real(8)   o2_w_root          ! Dry weight per root length (kg/m)
      real(8)   o2_w_root_z0       ! Root weight at depth
      real(8)   o2_soil_temp       ! Soil temperature (K)
      real(8)   o2_sat_water_cont  ! Saturated water content
      real(8)   o2_gas_filled_porosity ! Gas-filled porosity
      real(8)   o2_d_o2inwater     ! O2 diffusion in water
      real(8)   o2_d_root          ! Diffusion in root
      real(8)   o2_d_soil          ! Soil diffusion
      real(8)   o2_perc_org_mat    ! Organic matter percentage
      real(8)   o2_soil_density    ! Soil density (kg/m3)
      real(8)   o2_depth           ! Compartment thickness (m)
      real(8)   o2_shape_factor_microbialr ! Shape factor microbial resp
      real(8)   o2_root_radius     ! Root radius (m)
      real(8)   o2_r_microbial_z0  ! Microbial respiration rate
      real(8)   o2_waterfilm_thickness ! Water film thickness
      real(8)   o2_bunsencoeff     ! Bunsen coefficient
      real(8)   o2_c_min_micro     ! Min O2 for microbial resp
      real(8)   o2_c_macro         ! Macropore O2 conc
      real(8)   o2_ctopnode        ! Top node O2 conc
      logical   o2_initialized     ! Initialization flag
      
      ! [SS-BMI2 Task 5] retired 2026-05-13 — moved to state%timecontrol (ADR 0041)
      ! integer   period             ! moved to state%timecontrol%period
      ! integer   swheader           ! moved to state%timecontrol%swheader
      ! integer   swodat             ! moved to state%timecontrol%swodat
      ! integer   swres              ! moved to state%timecontrol%swres
      ! integer   swscre             ! moved to state%timecontrol%swscre
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%dt (ADR 0041)
      ! real(8)   dt                 ! moved to state%timecontrol%dt
      ! [SS-BMI2 Task 5] retired 2026-05-13 — moved to state%timecontrol (ADR 0041)
      ! real(8)   dtmax              ! moved to state%timecontrol%dtmax
      ! real(8)   dtmin              ! moved to state%timecontrol%dtmin
      real(8)   outdat(maout)      ! Array with output dates for water and solute balances
      real(8)   outdatint(maout)   ! Array with intermediate output dates
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%* (ADR 0041)
      ! real(8)   outper             ! moved to state%timecontrol%outper
      ! real(8)   t                  ! moved to state%timecontrol%t
      ! real(8)   t1900              ! moved to state%timecontrol%t1900
      ! real(8)   tcum               ! moved to state%timecontrol%tcum
      ! [SS-BMI2 Task 5] retired 2026-05-13 — moved to state%timecontrol (ADR 0041)
      ! real(8)   tend               ! moved to state%timecontrol%tend
      ! real(8)   tstart             ! moved to state%timecontrol%tstart
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%* (ADR 0041)
      ! logical   flbaloutput        ! moved to state%timecontrol%flbaloutput
      ! logical   fldayend           ! moved to state%timecontrol%flDayEnd
      ! logical   fldaystart         ! moved to state%timecontrol%flDayStart
      ! logical   fldecdtmin         ! moved to state%timecontrol%fldecdtmin
      ! logical   fldtmin            ! moved to state%timecontrol%fldtmin
      ! logical   fldtreduce         ! moved to state%timecontrol%fldtreduce
      ! logical   flheader           ! moved to state%timecontrol%flheader
      ! logical   floutput           ! moved to state%timecontrol%floutput
      ! logical   flrunend           ! moved to state%timecontrol%flRunEnd
      logical   flSwapShared       ! Flag to indicate the shared simultaneous simulation with other applications
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flYearStart (ADR 0041)
      ! logical   flyearstart        ! moved to state%timecontrol%flYearStart
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%date (ADR 0041)
      ! character(len=11) date       ! moved to state%timecontrol%date
      character(len=16) outfil     ! Name of output file
      character(len=80) pathwork   ! Path to work directory
      character(len=80) project    ! Name of project

! --- meteo variables
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%daymeteo (ADR 0041)
      ! integer   daymeteo           ! moved to state%timecontrol%daymeteo
      integer   daynrfirst         ! First calendar day number for which meteorological data is available in current year
      integer   daynrlast          ! Last calendar day number for which meteorological data is available in current year
      integer   detrecord(nmetfile) ! Record number of meteo file with detailed meteo data (-)
      integer   irectotal          ! Total record number with detailed meteo input for new weather file (-)
      integer   nmetdetail         ! Number of detailed records for ET and rainfall per day (-)
      integer   nmrain             ! Number of rain event records (-)
      integer   nofd               ! number of days for running average Tmin (-)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%rainrec (ADR 0041)
      ! integer   rainrec            ! moved to state%timecontrol%rainrec
      integer   swdivide           ! Switch on division ET into E and T: 0 = according to the SWAP traditional way; 1 = according to PMdirect
      integer   swetr              ! Switch: 0 = use daily meteorological basic data; 1 = use daily Etref values
      integer   swetsine           ! Switch: 0 = Tp and Ep uniform during a day; 1 = Tp and Ep are distributed as sine waves during a day
      integer   swinter            ! Switch for interception method: 0 = no interception; 1 = agricultural crops; 2 = trees and forests
      integer   swmetdetail        ! Switch: 0 = daily meteorological records; 1 = detailed records for both ET and rainfall
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%swmeteo (ADR 0041)
      ! integer   swmeteo            ! moved to state%timecontrol%swmeteo
      integer   swrain             ! Switch: 0 = use daily rain amounts; 1 = use daily amounts + mean intensity; 
                                   !         2 = use daily amounts + duration; 3 = use detailed rainfall data from separate file
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%* (ADR 0041)
      ! integer   wrecord            ! moved to state%timecontrol%wrecord
      ! integer   yearmeteo          ! moved to state%timecontrol%yearmeteo
      integer   ad(mrain)          ! Array with day numbers in meteo file
      integer   am(mrain)          ! Array with month numbers in meteo file
      real(8)   aetr(366)          ! Array with daily ETref input data (L/T)
      real(8)   ahum(366)          ! Array with daily humidity input data (M/L/T2)  
      ! ========================================================================
      ! [SS-ATM] retired 2026-05-11 — atmosphere subsystem migrated to state%atmosphere
      !   See: ADR 0037 (A-2.7)
      ! ========================================================================
      ! real(8) :: aintcdt   ! Interception flux of ONLY Rain during iteration timesteps (L/T)
      real(8)   alt                ! Altitude of meteorological station (L)
      real(8)   altw               ! Height of wind speed measurement (L)
      real(8)   angstroma          ! first  angstrom coefficient [-]
      real(8)   angstromb          ! second angstrom coefficient [-]
      real(8)   arad(366)          ! Array with daily radiation input data (M/T2)
      real(8)   arai(366)          ! Array with daily precipitation sum input data (L/T)
      real(8)   atav(96)           ! In case of detailed weather input, air temperature of each weather record (L/T)
      ! real(8) :: atmdem    ! Atmospheric demand = daily potential transpiration of a dry crop (L/T)
      real(8)   atmin7(7)          ! Array with minimum temperatures of last week (C)
      real(8)   atmn(366)          ! Array with daily minimum temperature input data (  )
      real(8)   atmx(366)          ! Array with daily maximum temperature input data (  )
      real(8)   awin(366)          ! Array with daily wind speed input data (L/T)
      ! real(8) :: caintc    ! Cumulative amount of rainfall interception (L)
      ! real(8) :: cevap     ! Cumulative amount of actual soil evaporation (L)
      ! real(8) :: cgrai     ! Cumulative amount of gross precipitation (L)
      ! real(8) :: cnrai     ! Cumulative amount of net precipitation (L)
      ! real(8) :: cpeva     ! Cumulative amount of potential soil evaporation (L)
      ! real(8) :: cptra     ! Cumulative amount of potential transpiration (L)
      real(8)   daylp              ! Photoperiodic daylength in hours (T)
      real(8)   dethum(nmetfile)   ! Array with detailed humidity input data (M/L/T2) 
      real(8)   detrad(nmetfile)   ! Array with detailed radiation input data (M/T2
      real(8)   detrain(nmetfile)  ! Array with detailed precipitation sum input data (L/T)
      real(8)   dettav(nmetfile)   ! Array with detailed temperature input data (  )
      real(8)   dettime(nmetfile)  ! Array with dates of detailed meteo input data
      real(8)   detwind(nmetfile)  ! Array with detailed wind speed input data (L/T)   
      real(8)   dtEventRain        ! Time step length for next precipitation event (T)
      ! real(8) :: empreva   ! Reduced soil evaporation flux according to empirical functions (L/T)
      real(8)   epot(96)           ! In case of detailed weather input, calculated Epot of each weather record (L/T)
      real(8)   cfevappond         ! Parameter equal to the ratio ponding layer evaporation / ETref (-)
      real(8)   finterception      ! Ratio net / gross rain flux in case of detailed rainfall data (-)
      ! real(8) :: fprecnosnow  ! Ratio rain (excl. snow and rain on snow) / gross rain flux (-)
      ! real(8) :: grai        ! Daily gross rain flux (L/T), without rain on snow
      ! real(8) :: graidt      ! Gross precipitation flux during iteration timesteps (L/T)
      real(8)   grain(96)          ! In case of detailed weather input, gross rain flux of each weather record (L/T)
      ! real(8) :: ievap       ! Intermediate amount of actual soil evaporation (L)
      ! real(8) :: inrai       ! Intermediate amount of net precipitation (L)
      ! real(8) :: ipeva       ! Intermediate amount of potential soil evaporation (L)
      ! real(8) :: iptra       ! Intermediate amount of potential transpiration (L)
      real(8)   lat                ! Latitude of meteorological station (degrees)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%metperiod (ADR 0041)
      ! real(8)   metperiod          ! moved to state%timecontrol%metperiod
      ! real(8) :: nraida      ! Daily average net precipitation flux (L/T)
      ! real(8) :: nraidt      ! Net precipitation flux during iteration timesteps (L/T)
      real(8)   nrain(96)          ! In case of detailed weather input, calculated netto rain of each weather record (L/T)
      ! real(8) :: peva        ! Potential soil evaporation flux (L/T)
      ! real(8) :: pevaday     ! Potential soil evaporation of one day (L)
      ! real(8) :: ptra        ! Potential transpiration flux (L/T)
      ! real(8) :: ptraday     ! Potential transpiration of one day (L)
      real(8)   rad                ! Global solar radiation (J/m2/d)
      real(8)   rainamount(mrain)  ! Array with short duration rainfall amounts (L)
      real(8)   raintab(60)        ! Array with mean rainfall intensity (L/T) as function of time (T)
      real(8)   rainfluxarray(mrain) ! Array with short duration rainfall intensities (L/T)
      real(8)   raintimearray(mrain) ! Array with times (T) at which rainfall intensity changes
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%rh (GR-ATM arc)
      ! real(8)   rh                 ! Relative air humidity (-)
      real(8)   tav                ! Average air temperature of a day (oC)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%tavd (GR-ATM arc)
      ! real(8)   tavd               ! Average air temperature during day time (oC)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%timjan1 (ADR 0041)
      ! real(8)   timjan1            ! moved to state%timecontrol%timjan1
      real(8)   tmn                ! Minimum air temperature of current day (oC)
      real(8)   tmnr               ! Average of minimum air temperature during past 7 days (oC)
      real(8)   tmx                ! Maximum air temperature of current day (oC)
      real(8)   tpot(96)           ! In case of detailed weather input, calculated Tpot of each weather record (L/T)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%tra (ADR 0038)
      ! real(8)   tra                ! Actual transpiration flux (L/T)
      real(8)   wet(366)           ! Fraction of each day the crop is wet (L)    
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%* (ADR 0041)
      ! logical   fletsine           ! moved to state%timecontrol%fletsine
      ! logical   flmeteodt          ! moved to state%timecontrol%flmeteodt
      ! logical   flmetdetail        ! moved to state%timecontrol%flmetdetail
      ! logical   flrainintens       ! moved to state%timecontrol%flrainintens
      ! logical   flupdmetdet        ! moved to state%timecontrol%flUpdMetDet
      character(len=200) metfil    ! Name of meteorological input file
      character(len=80) pathatm    ! Path to folder with meteorological input files
      ! CSV meteo cache: pre-loaded by adapter, sliced per year by MeteoCSVYear.
      ! Column layout (daily): 1=date, 2=rad(kJ/m2/d), 3=tmin, 4=tmax, 5=hum, 6=wind, 7=rain, 8=etref, 9=wet
      integer :: nmetcsv = 0
      real(8), dimension(:,:), allocatable :: metcsv_dat
      ! Detail CSV cache (swmetdetail=1): 7 columns per ADR 0014.
      ! 1=datetime(frac days since JD1900), 2=record, 3=rad(kJ/m2/d),
      ! 4=temp(C), 5=hum(kPa), 6=wind(m/s), 7=rain(mm)
      integer :: nmetcsv_det = 0
      real(8), dimension(:,:), allocatable :: metcsv_det
      ! Rain events CSV cache: pre-loaded by adapter, sliced per year by ReadRainEvents.
      ! Column layout: 1=datetime (fractional days since JD2415020), 2=amount (mm)
      integer :: nraincsv = 0
      real(8), dimension(:,:), allocatable :: raincsv_dat
!   - atmosphere SAVE variable state (refactored from local SAVE)
      real(8)   tsunrise_atm       ! Time of sunrise (fraction of day) - from meteodt.f90 ETSine
      real(8)   tsunset_atm        ! Time of sunset (fraction of day) - from meteodt.f90 ETSine  
      integer   nod10_cn           ! Node at -10cm for CN runoff method - from meteoday.f90 CNmethod
      integer   icn_atm            ! Current position in CN time table - from meteoday.f90 CNmethod
      real(8)   z10_cn             ! Depth to node 10 for CN method - from meteoday.f90 CNmethod
!   - meteo output variables for PEARL
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_etr (GR-ATM arc)
      ! real(4)   out_etr            ! Reference evapotranspiration  of current day (m/d)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_hum (GR-ATM arc)
      ! real(4)   out_hum            ! Air humidity  of current day (kPa)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_rad (GR-ATM arc)
      ! real(4)   out_rad            ! Global solar radiation (KJ/m2)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_tmn (GR-ATM arc)
      ! real(4)   out_tmn            ! Minimum air temperature of current day (oC)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_tmx (GR-ATM arc)
      ! real(4)   out_tmx            ! Maximum air temperature of current day (oC)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_wet (GR-ATM arc)
      ! real(4)   out_wet            ! Rainfall duration of current day (d)
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%out_win (GR-ATM arc)
      ! real(4)   out_win            ! Average wind speed of current day (m/s)

! --- irrigation variables
      integer   swirg              ! Switch for output file *.IRG with irrigation output: 0 = no; 1 = yes
      integer   irg                ! Internal number of irrigation output file *.IRG
      integer   irrigevent         ! Switch: 0 = no irrigation; 1 = fixed irrigation event; 2 = scheduled irrigation event
      integer   irtype(mairg)      ! Type of fixed irrigation: 0 = sprinkling irrigation; 1 = surface irrigation
      integer   isua               ! Switch for type of irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer   isuas              ! Switch for type of scheduled irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer   nirri              ! Number of irrigation event
      integer   phormc             ! Swith for 5th irrigation criterion: 0 = use pressure head; 1 = use water content
      integer   schedule           ! Switch for simulation of irrigation scheduling: 0 = no, 1 = yes
      integer   swirfix            ! Switch for fixed irrigation: 0 = no applications prescribed; 1 = applications are prescribed
      integer   swcirrthres        ! Switch to allow over irrigation when a conc-threshold is exceeded: 0 = no; 1 = yes/allowed
      real(8)   cirrs              ! Solute concentration of irrigation water (M/L3)
      real(8)   cirrthres          ! Threshold value (M/L3) indicating the concentration that initiates over irrigation
      real(8)   dcrit              ! Depth (L) of sensor for soil water pressure head or water content
      real(8)   ditab(14)          ! Array with amount of under- or over-irrigation (L) as function of crop development stage
      real(8)   dwatab(14)         ! Array with maximum amounts of water depleted as function of crop development stage
      real(8)   fidtab(14)         ! Array with prescribed fixed irrigation depth (L) as function of crop development stage
      real(8)   gird               ! Gross irrigation depth (L)
      real(8)   hcritab(14)        ! Array with minimum soil water pressure heads (L) as function of crop development stage
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   igird              ! Intermediate depth of gross irrigation (L)
      ! real(8)   inird              ! Intermediate depth of net irrigation (L)
      real(8)   irconc(mairg)      ! Array with irrigation concentrations (M/L3) in case of fixed irrigation
      real(8)   irdate(mairg)      ! Array with fixed irrigation dates
      real(8)   irdepth(mairg)     ! Array with fixed irrigation depths (L)
      real(8)   nird               ! Net irrigation depth (L)
      real(8)   perirrsurp         ! percentage (-) of the scheduled irrigation depths that may be over irrigated
      real(8)   raithreshold       ! Threshold value (L) indicating the amount of rainfall which is substracted from scheduled irrigation depths
      real(8)   rawtab(14)         ! Array with minimum of readily available water as function of crop development stage
      real(8)   tawtab(14)         ! Array with minimum of totally available water as function of crop development stage
      real(8)   tcritab(14)        ! Array with minimum volumetric soil water contents as function of crop development stage
      real(8)   tstairrig          ! Date after which scheduled irrigation is allowed
      real(8)   tendirrig          ! Date after which scheduled irrigation is NOT allowed
      real(8)   treltab(14)        ! Array with minimum of ratio actual/potential transpiration as function of crop development stage
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%* (ADR 0041)
      ! logical   flheadirg          ! moved to state%timecontrol%flheadirg
      ! logical   flirrigate         ! moved to state%timecontrol%flIrrigate
      ! logical   flIrg1Start        ! moved to state%timecontrol%flIrg1Start
      logical   FlIrrigationOutput ! Flag indication irrigation output
      integer   dayfix             ! days since last irrigation event

! --- start of growth grass      
      integer   swtsum             ! start of growth grass,1=TSUM200;2=TSOIL
      integer   tsumtime           ! time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      real(8)   tsumtemp           ! temperature limit to initiate grass growth  [0.0..20.0 grC, R]
      real(8)   tsumdepth          ! depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
      
! --- crop variables
      integer   swCrop             ! Switch for simulating crop: 0 = No; 1 = Yes
      logical   flCropCalendar     ! Flag indicating that crop season is active (but currently might be bare or cropped)
      logical   flCropEmergence    ! Flag indicating period from crop emergence until harvest
      logical   flCropHarvest      ! Flag indicating period from crop harvest until the end of crop season
      logical   flCropReadFile     ! Flag indicating reading of input.crp
      logical   flCropOpenFile     ! Flag indicating to create output.crp
      logical   flCropOutput       ! Flag indicating writing of output.crp

      integer   croptype(macrop)   ! Switch for type of crop model: 1 = simple; 2 = general detailed model; 3 = detailed for grass
      integer   swcrp              ! Switch for output file *.CRP with daily crop output: 0 = no; 1 = yes
      integer   crp                ! Internal number of crop output file *.CRP
      integer   daycrop            ! Number of days that a crop exists
      integer   icrop              ! Current crop number
      integer   daygrowth          ! Number of days that grass is growing after management event or emergence, actual run (d)
      integer   daygrowthpot       ! Number of days that grass is growing after management event or emergence, potential run (d)
      integer   idaysgraz          ! Day number of grass grazing, actual run (d)
      integer   idaysgrazpot       ! Day number of grass grazing, potential run (d)
      integer   idev               ! Switch for length of growth period in case of simple crop: 1 = fixed; 2 = depends on temperature sum
      integer   idregr             ! Number of days for regrowth of grassland, actual growth
      integer   idregrpot          ! Number of days for regrowth of grassland, potential growth
      integer   idsl               ! Switch for crop development before anthesis: 0 = depends on temperature; 
                                   !   1 = depends on temperature and day length; 2 = depends on temperature, day length and vernalisation factor
      integer   iharvest           ! Grass harvest number of grass crop when harvest dates are fixed
      integer   ilvold             ! Age of oldest leaf (d) of actual crop
      integer   ilvoldpot          ! Age of oldest leaf (d) of potential crop
      integer   iseqgm             ! Counter in sequence of grass grazing and mowing, actual crop
      integer   iseqgmpot          ! Counter in sequence of grass grazing and mowing, potential crop
      integer   noddrz             ! Compartment number at bottom root zone (-)
      integer   seqgrazmow(366)    ! Sequence of grass grazing and mowing, actual crop
      integer   seqgrazmowpot(366) ! Sequence of grass grazing and mowing, potential crop
      integer   swcf               ! Switch for simple crop: 1 = crop factor is input; 2 = crop height is input
      integer   swdrought          ! Switch for drought stress (1 = Feddes et al., 1978; 2 = De Jong van Lier et al., 2008)
      integer   swgc               ! Switch for simple crop: 1 = leaf area index is input; 2 = soil cover fraction is input
      integer   swjarvis           ! Switch for compensation of root water uptake stress according to Jarvis (1989) (deprecated from 4.1.73; replaced by swcompensate)
      integer   swcompensate       ! Switch for method of compensation of root water uptake stress
      real(8)   alphacrit          ! Critical stress index for compensation of root water uptake (-)
      real(8)   dcritrtz           ! Threshold for rootzone to start compensation of root water uptake; Walsum (cm)
      integer   swstressor         ! Switch for stressor to compensate (1 = all stressors (default); 2 = drought stress, 3 = oxygen stress, 4 = salinity stress; 5 = frost stress)
      integer   swrootradius       ! Switch for root radius to calculate Oxygen stress: 1 = calculate, 2 = input
      integer   swsalinity         ! Switch for salinity stress: 0 = no stress; 1 = Maas and Hoffman (1977); 2 = Osmotic head
      real(8)   atmtr              ! Daily atmospheric transmission (-)
      real(8)   agerm              ! Coefficient a  of germination
      real(8)   cgerm              ! Coefficient c  of germination
      real(8)   bgerm              ! Coefficient b  of germination
      real(8)   adcrh              ! Level of high atmospheric demand (L/T)
      real(8)   adcrl              ! Level of low atmospheric demand (L/T)
      real(8)   air_filled_root_por ! Air filled root porosity [0..1.0 -, R]
      real(8)   albedo             ! Crop reflection coefficient (-)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%alpJvLier (ADR 0036)
      ! real(8)   alpJvLier          ! Uniform drought reduction factor based on concept Jong van Lier (-)
      real(8)   amaxtb(30)         ! Maximum CO2 assimilation rate (kg/ha/hr) as function of development stage (-)
      real(8)   pgasspot           ! Assimilation rate after nitrogen stress and maximum attainable yield, potential crop growth
      real(8)   pgass              ! Assimilation rate after nitrogen stress and maximum attainable yield, actual crop growth
      integer   swpotrelmf         ! Calculation of potential yield
      real(8)   relmf              ! Management factor (attainable yield)
      real(8)   avevaptb(2*magrs)  ! Gash interception model: average evaporation intensity during shower (-) as function of time (T)
      real(8)   avprectb(2*magrs)  ! Gash interception model: average rainfall intensity (-) as function of time (T)
      real(8)   cf                 ! Crop factor (-)
      real(8)   cfeic              ! Crop factor wet (-)
      real(8)   gctb(2*magrs)      ! Array with either LAI (m2/m2) or Soil Cover Fraction (m2/m2) 
      real(8)   cftb(2*magrs)      ! Array with either crop factors (-) or crop height (L) as function of development stage
      real(8)   cfeictb(2*magrs)   ! Array with crop factors, wet (-) 
      real(8)   ch                 ! Crop height (cm)
      real(8)   chtb(2*magrs)      ! Array with crop heights (cm) as function of development stage
      real(8)   c_mroot            ! Maintenance coefficient of root [0.0..1.0 kg O2/kg/d, R]
      real(8)   cofab              ! Interception coefficient Von Hoyningen-Hune and Braden (L)
      real(8)   cropend(macrop)    ! Array with crop end dates
      real(8)   cropstart(macrop)  ! Array with crop start dates
      real(8)   cumdens(202)       ! Cumulative root density as function of relative soil depth (-)
      real(8)   cuptgraz           ! Cumulative dry weight of grass consumed with animal grazing for actual run (kg/ha)
      real(8)   cuptgrazpot        ! Cumulative dry weight of grass consumed with animal grazing for potential run (kg/ha)
      real(8)   cvl                ! Efficiency of assimilate conversion into leaves (kg/kg)
      real(8)   cvo                ! Efficiency of assimilate conversion into storage organs (kg/kg)
      real(8)   cvr                ! Efficiency of assimilate conversion into roots (kg/kg)
      real(8)   cvs                ! Efficiency of assimilate conversion into stems (kg/kg)
      real(8)   cwdm               ! Dry weight of dead and living plant organs (kg/ha)
      real(8)   cwdmpot            ! Dry weight of dead and living plant organs for potential growth (kg/ha)
      real(8)   difpp              ! Diffuse irradiation perpendicular to direction of light (J/m2/s)
      real(8)   dlc                ! Shortest day length (T) for any crop development
      real(8)   dlo                ! Minimum day length (T) for optimal crop development
      real(8)   dry_mat_cont_roots ! Dry matter content of roots [0..1.0 -, R]
      real(8)   dsinbe             ! Daily total of effective solar height (s)
      real(8)   dtsmtb(30)         ! Increase of temperature sum (oC) as function of daily average temperature (oC)
      real(8)   dvs                ! Crop development stage (-)
      real(8)   dvsend             ! Crop development stage at harvest
      real(8)   dwlv               ! Dry weight of plant leafs of actual crop (kg/ha)
      real(8)   dwlvCrop           ! Dry weight of deceased leafs that remain attached to plant (kg/ha)
      real(8)   dwlvSoil           ! Dry weight of deceased leafs allocated to soil (kg/ha)
      real(8)   dwlvpot            ! Dry weight of plant leafs of potential crop (kg/ha)
      real(8)   dwrt               ! Dry weight of plant roots of actual crop (kg/ha)
      real(8)   dwrtpot            ! Dry weight of plant roots of potential crop (kg/ha)
      real(8)   dwso               ! Dry weight of plant storage organs of actual crop (kg/ha)
      real(8)   dwst               ! Dry weight of plant stem of actual crop (kg/ha)
      real(8)   dwstpot            ! Dry weight of plant stem of potential crop (kg/ha)
      real(8)   eff                ! Light use efficiency of a leaf (kg CO2 / J adsorbed)
      real(8)   f_senes            ! Reduction factor for senescence, used for maintenance respiration [0..1.0 -, R]
      real(8)   fltb(30)           ! Fraction of total dry matter increase partitioned to the leaves (-) as function of dvs
      real(8)   fotb(30)           ! Fraction of total dry matter increase partitioned to the storage organs (-) as function of dvs
      real(8)   frtb(30)           ! Fraction of total dry matter increase partitioned to the roots (-) as function of dvs
      real(8)   fstb(30)           ! Fraction of total dry matter increase partitioned to the stems (-) as function of dvs
      real(8)   gasst              ! Total gross assimilation for actual crop (kg/ha)
      real(8)   gasstpot           ! Total gross assimilation for potential crop (kg/ha)
      real(8)   gc                 ! Ground cover in case of a crop
      real(8)   HarLosOrm_tot      ! Harvest losses added to the soil (roots+fr_shoots_fr_stems+fr_stor.organs) at harvest (kg/ha DM) 
      real(8)   hdrygerm           ! Criterium Hdry of germination
      real(8)   hwetgerm           ! Criterium Hwet of germination
      real(8)   hlim1              ! Pressure head above which root water uptake stops (L)
      real(8)   hlim2l             ! Pressure head below which optimum water uptake starts for sub layer (L)
      real(8)   hlim2u             ! Pressure head below which optimum water uptake starts for top layer (L)
      real(8)   hlim3h             ! Pressure head below which water uptake reduction starts at high Tpot (L)
      real(8)   hlim3l             ! Pressure head below which water uptake reduction starts at low Tpot (L)
      real(8)   hlim4              ! Wilting point, no root water uptake at lower soil water pressure heads (L)
      real(8)   kdif               ! Extinction coefficient for diffuse visible light (-)
      real(8)   kdir               ! Extinction coefficient for direct visible light (-)
      real(8)   lai                ! Leaf area index
      real(8)   laiem              ! Leaf area index (-) at crop emergence
      real(8)   laiexp             ! Leaf area index (-) in exponential growth stage of actual crop
      real(8)   laiexppot          ! Leaf area index (-) in exponential growth stage of potential crop
      real(8)   glaiex             ! increase in leaf area during exponential growth of actual crop
      real(8)   glaiexpot          ! increase in leaf area during exponential growth of potential crop
      real(8)   laimax             ! Maximum leaf area index reached during growth of actual grass crop (-)
      real(8)   laipot             ! Leaf area index for potential run (-)
      real(8)   lv(366)            ! Array with leaf weight (kg/ha) as function of crop day number of actual crop
      real(8)   lvpot(366)         ! Array with leaf weight (kg/ha) as function of crop day number of potential crop
      real(8)   lvage(366)         ! Array with leaf age (d) as function of crop day number of actual crop
      real(8)   lvagepot(366)      ! Array with leaf age (d) as function of crop day number of potential crop
      real(8)   max_resp_factor    ! Ratio root total respiration / maintenance respiration [1..5.0 -, R]
      real(8)   mowrest            ! Dry weight of above ground grass (leaves + stems) after mowing (kg/ha)
      real(8)   dewrest            ! Dry weight of above ground grass (leaves + stems) after dewooling (kg/ha)
      real(8)   mrest              ! Total maintenance respiration for actual crop (kg/ha)
      real(8)   mrestpot           ! Total maintenance respiration for potential crop (kg/ha)
      real(8)   mrftb(2*magrs)     ! Array with ratio root total respiration / maintenance respiration as function of DVS (kg/m3)
      real(8)   perdl              ! Maximum relative death rate of leaves due to water stress (/T)
      real(8)   pfreetb(2*magrs)   ! Gash interception model: free throughfall coefficient (-) as function of time (T)
      real(8)   pstemtb(2*magrs)   ! Gash interception model: stem flow coefficient (-) as function of time (T)
      real(8)   siccaptb(2*magrs)  ! NHI interception model: interception capacity as a function of time (T)
      real(8)   fimin              ! start-up saturation fraction for relative interception evaporation (-)
      real(8)   siccapact          ! interceprion storage capacity of canopy (cm)
      ! real(8) :: sicact        ! amount of water stored on canopy (cm) — [SS-ATM] retired 2026-05-11
      real(8)   siccaplai          ! interception storage per unit of LAI (cm/LAI)
      real(8)   q10                ! Relative increase of respiration rate with temperature (/10 oC)
      real(8)   q10_microbial      ! Relative increase in microbial respiration at temperature increase of 10 �C [1.0..4.0 -, R]
      real(8)   q10_root           ! Relative increase in root respiration at temperature increase of 10 �C [1.0..4.0 -, R]
      real(8)   rdrrtb(30)         ! Array with relative death rates of roots (kg/kg/d) as function of development stage (-)
      real(8)   rdrstb(30)         ! Array with relative death rates of stems (kg/kg/d) as function of development stage (-)
      real(8)   reltr              ! relative transpiration factor that reduces crop growth (-)
      real(8)   rfsetb(30)         ! Reduction factor for senescence (-) as function of development stage (-)
      real(8)   rgrlai             ! Maximum relative increase in leaf area index (/T)
      real(8)   rid                ! Real day number of detailed grass crop (d)
      real(8)   rml                ! Relative maintenance respiration rate of leaves (kg CH2O)/kg/d)
      real(8)   rmo                ! Relative maintenance respiration rate of storage organs (kg CH2O)/kg/d)
      real(8)   rmr                ! Relative maintenance respiration rate of roots (kg CH2O)/kg/d)
      real(8)   rms                ! Relative maintenance respiration rate of stems (kg CH2O)/kg/d)
      real(8)   rootcoefa          ! Defines relative distance at which mean soil water content occurs between roots
      real(8)   rooteff            ! Root system efficiency factor [0..1.0 -, R]
      real(8)   rootradius         ! Root radius drought stress (cm)
      real(8)   root_radiusO2      ! Root radius oxygen stress (m)
      real(8)   rsc                ! Minimum canopy resistance of dry crop (T/L)
      real(8)   rsw                ! Canopy resistance of intercepted water (T/L)
      real(8)   scanopytb(2*magrs) ! Gash interception model: storage capacity of canopy (-) as function of time (T)
      real(8)   shape_factor_rootr ! Shape factor for exponential decrease of root respiration rate with depth [0..1.0 -, R]
      real(8)   sla(366)           ! Array with specific leaf area (ha/kg) as function of crop day number of actual crop
      real(8)   slapot(366)        ! Array with specific leaf area (ha/kg) as function of crop day number of potential crop
      real(8)   slatb(30)          ! Array with specific leaf area (ha/kg) as function of development stage
      real(8)   spa                ! Specific pod area (ha/kg)
      real(8)   span               ! Life span of leaves at optimum conditions (T)
      real(8)   spec_weight_root_tissue ! Specific weight of non-airfilled root tissue [0.d0..1.d5 kg root/m3 root, R]
      real(8)   specific_resp_humus ! Respiration rate of humus at 25 �C [0.0..1.0 kg O2/kg C/d, R] 
      real(8)   srl                ! Specific root length [0.d0..1.d10 m root/kg root, R]      
      real(8)   ssa                ! Specific stem area (ha/kg)
      real(8)   tadw               ! Dry weight of plant minus roots of actual growth (kg/ha)
      real(8)   tadwpot            ! Dry weight of plant minus roots of potential growth (kg/ha)
      real(8)   tagp               ! Dry weight of dead and living grass organs (kg/ha)
      real(8)   tagppot            ! Dry weight of dead and living grass organs for potential run (kg/ha)
      real(8)   tagpt              ! Dry weight of harvested grass (kg/ha)
      real(8)   tagptpot           ! Dry weight of harvested grass for potential run (kg/ha)
      real(8)   tbase              ! Lower threshold temperature for ageing of leaves (oC)
      real(8)   tdwi               ! Initial total crop dry weight (kg/ha)
      real(8)   tmnftb(30)         ! Reduction factor for maximum assimilation rate (-) as function of minimum day temperature (oC)
      real(8)   tmpftb(30)         ! Reduction factor for maximum assimilation rate (-) as function of average day temperature (oC)
      real(8)   tsum               ! Temperature sum from cropstart to cropend (oC)
      real(8)   tsumam             ! Temperature sum from anthesis to maturity (oC)
      real(8)   tsumea             ! Temperature sum from emergence to anthesis (oC)
      real(8)   tsumemeopt         ! Temperature sum for crop emergence under optimal conditions
      real(8)   tsumgerm           ! Temperature sum during germination
      real(8)   TBASEM             ! Lower threshold temp. for emergence (C)
      real(8)   TEFFMX             ! max. eff. temp. for emergence (C)

      real(8)   var_a              ! Variance of root radius [0.d0..1.d0 -, R]
      real(8)   w_root_ss          ! Dry weight of roots at soil surface [0.0..10.0 kg/m3, R]
      real(8)   wiltpoint          ! Minimum pressure head at interface soil-root (cm)

      real(8)   wlv                ! Dry weight of plant leaves (kg/ha)
      real(8)   wlvpot             ! Dry weight of plant leaves for potential growth (kg/ha)
      real(8)   wrtb(2*magrs)      ! Array with dry weight of root at soil surface as function of DVS (kg/m3)
      real(8)   wrt                ! Dry weight of plant root (kg/ha)
      real(8)   wrtpot             ! Dry weight of plant root for potential growth (kg/ha)
      real(8)   wrtmin             ! Minimum dry weight of plant root at relative depth (1% of the initial value)
      real(8)   gwrt               ! Growth of dry weight of plant root (kg/ha)
      real(8)   wso                ! Dry weight of storage organ (kg/ha)
      real(8)   wsopot             ! Dry weight of storage organ for potential growth (kg/ha)
      real(8)   wst                ! Dry weight of plant stem (kg/ha)
      real(8)   wstpot                ! Dry weight of plant stem for potential growth (kg/ha)
      logical   flanthesis         ! Flag indicating anthesis stage of a crop
      logical   flHarvest           ! Flag indicating that grass should be harvested, actual crop
      logical   flHarvestDay        ! Flag indicating that current day is harvest day
      logical   flHarvestpot        ! Flag indicating that grass should be harvested, potential crop
      logical   flGrazing           ! Flag indicating that cattle grazes the grass, actual run
      logical   flGrazingpot        ! Flag indicating that cattle grazes the grass, potential run
      character(len=40) cropfil(macrop)   ! Array with names of crop files
      character(len=80) pathcrop          ! Path to folder with crop input files
      character(len=200) inifil           ! Name of file with output data of previous day which is used for initialization

!     Harvest Grass
      real(8)   dmmowtb(20)             ! Array with threshold of mowing event as function of above ground dry matter
      real(8)   dmgrztb(20)             ! Array with threshold of grazing event as function of above ground dry matter
      real(8)   dateharvest(999)        ! Array with dates of mowing/grazing event
      real(8)   DelayRegrowthTab(2*100) ! Array with delay of regrowth as function of dry matter harvest
      real(8)   lsda(366)               ! Array with Lifestock density at grazing event
      real(8)   DaysGrazingtab(2*100)   ! Array with days of grazing as function lifestock density
      real(8)   UptGrazingtab(2*100)    ! Array with grazing uptake as function lifestock density
      real(8)   LossGrazingtab(2*100)   ! Array with grazing losses as function lifestock density
      real(8)   lossmowtab(2*100)       ! Array with extra dry matter losses during mowing event as function of pressure head
      real(8)   lossgrztab(2*100)       ! Array with extra dry matter losses during grazing event as function of pressure head
      
      real(8)   cropstartpot            ! Start of potential grass growth (reset after mowing or grazing event)
      real(8)   cropendpot              ! End of potential grass growth (reset after mowing or grazing event)
      real(8)   cropstartact            ! Start of actual grass growth (reset after mowing or grazing event)
      real(8)   cropendact              ! End of actual grass growth (reset after mowing or grazing event)
      logical   flhrvendpot             ! Flag indicating end of harvest event (potential)
      logical   flhrvendact             ! Flag indicating end of harvest event (actual)
      real(8)   pmowdm                  ! Total potential harvest by mowing at end of harvest event (kg/ha)
      real(8)   mowdm                   ! Total actual harvest by mowing at end of harvest event (kg/ha)
      real(8)   pgrzdm                  ! Total potential harvest by grazing at end of harvest event (kg/ha)
      real(8)   grzdm                   ! Total actual harvest by grazing at end of harvest event (kg/ha)
      real(8)   plossdm                 ! Total loss of potential harvest due to insufficient pressure head (kg/ha)
      real(8)   lossdm                  ! Total loss of actual harvest due to insufficient pressure head (kg/ha)
      
! --- rooting      
      integer   swrdc              ! Switch for calculation of root density (0: static; 1: dynamic)
      real(8)   rdctb(22)          ! Array with relative root density (-) as function of relative root depth (-)
      integer   swrd               ! Switch for development of root extension
      real(8)   rdtb(2*magrs)      ! root depth as function of development stage
      real(8)   rlwtb(22)          ! root depth as function of root biomass
      real(8)   wrtmax             ! maximum root weight
      integer   swdmi2rd           ! rooting depth limitation by relative dry matter increase (dmi/dmipot)
      real(8)   rdi                ! Initial rooting depth (L)
      real(8)   rri                ! Maximum daily increase of rooting depth (L/T)
      real(8)   rdc                ! Maximum rooting depth of particular crop (L)
      real(8)   rdmax              ! Maximum rooting depth in soil profile (L)
      real(8)   rd                 ! Rooting depth (L)
      real(8)   rdpot              ! Rooting depth for potential run (L)
      real(8)   rdm                ! Maximum rooting depth (minimum of soil profile and particular crop) (L)

! --- CO2
      logical   flCO2              ! Flag indicating correction of CO2
      real(8)   fco2amax           ! factor to correct AMAX for CO2
      real(8)   fco2eff            ! factor to correct EFF for CO2
      real(8)   fco2tra            ! factor to correct TRA for CO2
      real(8)   co2amaxtb(30)      ! table with factors to correct AMAX for CO2
      real(8)   co2efftb(30)       ! table with factors to correct EFF for CO2
      real(8)   co2tratb(30)       ! table with factors to correct TRA for CO2
      integer   co2year(mayrs)     ! table with years for which CO2 concentrations are given
      real(8)   co2ppm(mayrs)      ! table with CO2 concentrations (ppm), for each year in co2year

! --- vernalisation
      real(8)   verndvs            ! critical development stage after which the effect of vernalisation is halted [-]
      real(8)   vernsat            ! saturated vernalisation requirement [d]
      real(8)   vernbase           ! base vernalisation requirement [d]
      real(8)   vernrtb(30)        ! table with rate of vernalisation as function of tav [days/degrees]
!     only for bulb crops (tulips etc..)
      integer   swbulb             ! switch to enable simulation of bulb crops (-)
      real(8)   drbl               ! Death rate of actual bulb (kg/ha)
      real(8)   drblpot            ! Death rate of potential bulb (kg/ha)
      real(8)   dwbl               ! Dry weight of dead flowers (kg/ha)
      real(8)   dwblpot            ! Dry weight of dead flowers for potential growth (kg/ha)
      real(8)   fbl                ! Dry weight fraction partitioned to flowers (-)
      real(8)   fbltb(30)          ! Dry weight fractions partitioned to flowers (-)
      real(8)   pld
      real(8)   remoc
      real(8)   plwt               ! Dry weight of mother bulb (kg/ha)
      real(8)   plwti              ! Initial dry weight of mother bulb (kg/ha)
      real(8)   wbl                ! Dry weight of living flowers (kg/ha)
      real(8)   wblpot             ! Dry weight of living flowers for potential growth (kg/ha)

! --- Nitrogen: crop and soil management
      logical   flCropNut          ! Flag indicating simulation of crop nutrient stress
      logical :: flTillage = .false.   !! ADR 0020 call-site gate for DoTillage
      logical :: flSSDI    = .false.   !! ADR 0020 call-site gate for SSDI_irrigation
      real(8)   nmxlv(30)
      real(8)   nlue,anlv,anst,nmaxlv,nmaxst,nmaxrt
      real(8)   lrnr,lsnr,nni,rnflv,rnfst,frnx
      ! N-P-K nutrient parameters from cropwofost.nutrient (N1 of [nutrients] umbrella).
      ! Promoted from local-to-wofost-subroutine after spec brainstorming
      ! revealed they could not be reached from a config-load-time adapter
      ! while declared as locals. See ADR 0025.
      real(8)   nlai, nmaxso, npart, nfixf
      real(8)   nsla, rnfrt, tcnt
      real(8)   dvsnlt, dvsnt, rdrns, fntrt
      integer   ilnmxl
      real(8)   fraharlosorm_lv, fraharlosorm_st, fraharlosorm_so
      real(8)   fstr
      real(8)   amFERT             ! amount of applied Fertilizer (kg/ha/d N)

! --- tillage variables (bridge for tillage module)
      integer   till_swtill                       ! Switch: 0=no tillage, 1=tillage
      integer   till_i_n_model                    ! Switch for n-parameter treatment (1-3)
      integer   till_iRedist                      ! Redistribution type after MvG change
      integer   till_Ntill                        ! Number of tabulated tillage events
! ========================================================================
! [SS-TIL T-5] retired 2026-05-12 — tillage runtime-state globals moved to state%tillage
!   See: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
!        ADR 0039 (Task T-6)
! ========================================================================
!     integer   till_iTill                        ! Current tillage event index
      integer   till_Ntypes                       ! Number of tillage types
!     integer   till_MaxNumSoilHo                 ! Max soil horizons in tillage zone
!     integer   till_MaxNumSoilCP                 ! Max soil compartments in tillage zone
      real(8)   till_Max_Z_tillage                ! Max possible depth of tillage (cm)
      real(8), dimension(:), allocatable :: till_Date_tillage   ! Tillage dates
      real(8), dimension(:), allocatable :: till_Z_tillage      ! Tillage depths (cm)
      real(8), dimension(:), allocatable :: till_I_tillage      ! Tillage intensity (0-1)
      integer, dimension(:), allocatable :: till_Type_Tillage   ! Tillage type index
      integer, dimension(:), allocatable :: till_iType_Tillage  ! Tillage type identifier
      integer, dimension(:), allocatable :: till_iTT1           ! First position in type table
      integer, dimension(:), allocatable :: till_iTT2           ! Last position in type table
      real(8), dimension(:), allocatable :: till_TAB_Rho_tillage ! Bulk density after tillage
      real(8), dimension(:), allocatable :: till_TAB_Rho_cons   ! Consolidated bulk density
      real(8), dimension(:), allocatable :: till_TAB_K_R_cons   ! Consolidation rate constant
      real(8), dimension(:), allocatable :: till_TAB_Rho_match  ! Matching point density
      real(8), dimension(:), allocatable :: till_TAB_N_match    ! Matching point n-value
!     real(8), dimension(:), allocatable :: till_Rho_tillage    ! Post-tillage bulk density per layer
!     real(8), dimension(:), allocatable :: till_Rho_cons       ! Consolidated density per layer
!     real(8), dimension(:), allocatable :: till_Rho_last       ! Previous density per layer
!     real(8), dimension(:), allocatable :: till_K_R_cons       ! Consolidation rate per layer
!     real(8), dimension(:), allocatable :: till_Rho_match      ! Matching point density per layer
!     real(8), dimension(:), allocatable :: till_N_match        ! Matching point n per layer
!     real(8), dimension(:), allocatable :: till_Slope_match    ! Slope at matching point per layer
!     real(8)   till_sumDWC                       ! Sum of water content changes
!     real(8)   till_sumAvail1                    ! Available pore space
!     real(8)   till_sumAvail2                    ! Available water

! --- soilwater variables
      integer   iHWCKmodel(maho)   ! indicator what type of water retention and hydraulic conductivity model is used (per soil layer)
                                   ! 1 = MvG (default), 2 = exponential, 3 = MvG bi-modal
                                   ! 4-11: 8 versions of PDI model
                                   ! other types may be added in the future
      logical   BiModal(maho)      ! logical indicating whether chosen model is bi-modal or not
      logical   NoVap(maho)        ! logical indicating that NO vapour flow is to be considered in PDI K-model
      
      ! [SS-BMI2 Task 5] retired 2026-05-13 — moved to state%timecontrol (ADR 0041)
      ! integer   MaxIterTime        ! moved to state%timecontrol%MaxIterTime
      ! integer   MaxIt              ! moved to state%timecontrol%MaxIt
      integer   MaxBackTr
      integer   Itnumb(100,2)
      real(8)   CritDevh1Cp        ! Convergence criterium for Richards equation: relative difference in pressure heads (-)
      real(8)   CritDevh2Cp        ! Convergence criterium for Richards equation: absolute difference in pressure heads (L)
      real(8)   CritDevPondDt
      logical   fldumpconvcrit     ! flag to generate additional output about convergence-warnings from subr Headcalc
      logical   flwarn_hc          ! Headcalc warning flag (previously SAVE variable)
      integer   iwarn_hc           ! Headcalc warning counter (previously SAVE variable)
      integer   nstep_hc           ! Headcalc step counter (previously SAVE variable)
      integer   dev_cmb            ! Mass balance deviation file unit (previously SAVE in checkmassbal)
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%flksatexm (GR-BH arc)
      ! logical   flksatexm          ! flag Ksatexm variable present in input file
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%fluseksatexm (ADR 0038)
      ! logical   fluseksatexm(macp) ! flag per node: yes/no make use of Ksatexm (Ksat examined in lab or field) extension in h-range [-2,0]
      ! [SS-BMI2 Task 5] retired 2026-05-13 — moved to state%timecontrol (ADR 0041)
      ! logical   flMaxIterTime      ! moved to state%timecontrol%flMaxIterTime
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! logical   FlRunoff
      logical   swcaprise          ! flag to minimize cap.rise to rootzone (for experts only)
      logical   swcapriseoutput    ! flag to generate an output file with cap.rise to/form rootzone (for experts only)

      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%H0max/k1max/q0 (GR-BH arc)
      ! real(8)   h0max
      ! real(8)   k1max
      ! real(8)   q0

      integer   afo                ! Internal number of output file *.AFO with formatted hydrologic data for soil water quality models
      integer   aun                ! Internal number of output file *.AUN with unformatted hydrologic data for soil water quality models
      integer   bal                ! Internal number of output file *.BAL with overview of water balance
      integer   blc                ! Internal number of output file *.BLC with all water balance components in detail
      integer   bma                ! Internal number of output file *.BMA with detailed yearly water balance in case of macropores
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%bpegwl (ADR 0038)
      ! integer   bpegwl             ! Node at bottom of perched groundwater
      integer   botcom(maho)       ! Array with number of bottom compartments in each soil layer
      integer   dra                ! Internal number of drainage input file *.DRA
      integer   dramet             ! Switch for lateral drainage: 1 = table of flux - groundwater level; 2 = Hooghoudt or Ernst; 
                                   !                              3 = drainage/infiltration resistance
      integer   swinc              ! Switch for output file *.INC with daily incremental water balance data: 0 = no; 1 = yes
      integer   inc                ! Internal number of output file *.INC with incremental water balance data
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%indeks (ADR 0038)
      ! integer   indeks(macp)       ! Index denoting wetting or drying curve in case of hysteresis: 1 = wetting; -1 = drying
      integer   ipos               ! Switch for position of drain (see *.DRA input file for overview)
      integer   isoillay(maho)     ! Number of soil layer, starting with 1 at the soil surface
      ! [GR-BH Task 35] layer(macp) retired — moved to state%mesh%layer
      ! [SS-BMI2 Task 5] retired 2026-05-13 — moved to state%timecontrol%msteps (ADR 0041)
      ! integer   msteps             ! moved to state%timecontrol%msteps
      integer   ncomp(macp)        ! Array with number of compartments in each sublayer
      integer   nhead              ! Number of initial soil water pressure heads as provided in the input
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%nodgwl (ADR 0038)
      ! integer   nodgwl             ! Node directly above groundwater level
      integer   nod1lay(maho)      ! node nr of first node of each soil layer (from top to bottom)
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%nodfrostbot (ADR 0034)
      ! integer   nodfrostbot        ! Node nr of deepest node with frost conditions
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%npegwl (ADR 0038)
      ! integer   npegwl             ! Node directly above perched groundwater level
      ! [GR-BH Task 37] nrlevs retired — moved to state%drainage%nrlevs
      integer   nrstaring          ! Number of soil type [1..18] according to Staring series (Wosten et al., 2001)
      integer   nsublay            ! Number of sublayers in the soil profile
      integer   numbit             ! Iteration number for solving Richards equation
      integer   numlay             ! Number of (physical) soil layers
      ! [GR-BH Task 35] numnod retired — moved to state%mesh%numnod
      integer   numnodnew          ! Number of desired nodes for soil water quality models
      integer   numtab(macp)       ! Number of table entries of soil physical values for each model compartment
      integer   numtablay(maho)    ! Number of table entries of soil physical values for each soil layer
      integer   rot                ! Internal number of output file *.ROT with microscopic root water extraction data 
      integer   swstr              ! Switch for output file *.STR with daily stress factors output: 0 = no; 1 = yes
      integer   str                ! Internal number of output file *.STR with stress factors for transpiration
      integer   sw2                ! Switch for prescribed bottom flux: 1 = sine function; 2 = table
      integer   sw3                ! Switch for prescribed hydraulic head of deep aquifer: 1 = sine function; 2 = table
      integer   sw4                ! Switch for extra groundwater flux as function of time: 0 = no extra flux; 1 = include extra flux
      integer   swcsv              ! Switch for special CSV output specified by user; default = 0 ;if 1, requires InList_csv
      character(len=1024) InList_csv   ! character string with comma-separated list of variables for CSV output
      integer   swcsv_tz              ! Switch for special CSV output specified by user; default = 0 ;if 1, requires InList_csv
      character(len=1024) InList_csv_tz   ! character string with comma-separated list of variables for CSV output
      real(8)   tz_z1_z2(2)        ! Depth range for time-depth CSV output (default: top soil profile, bottom soil profile)
      integer   swafo              ! Switch for extra output file with formatted data for water quality models: 
                                   !       0 = no output; 1 = output to file *.AFO; 2 = output to file *.BFO
      integer   swaun              ! Switch for extra output file with unformatted data for water quality models: 
                                   !       0 = no output; 1 = output to file *.AUN; 2 = output to file *.BUN
      integer   swbal              ! Switch for output file with yearly water balance *.BAL: 0 = no; 1 = yes
      integer   swblc              ! Switch for output file with detailed yearly water balance *.BLC: 0 = no
      integer   swsba              ! Switch for output file with daily solute balance *.SBA: 0 = no; 1 = yes; 1 = yes
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%swbotb_runtime (GR-BH arc)
      ! integer   swbotb             ! Switch for bottom boundary condition (see *.SWP input file for overview)
      integer   swbotb3Impl        ! Switch for implicit solution with lower boundary option 3 (Cauchy): 0 = explicit, 1 = implicit
      integer   SwBotb3ResVert     ! Switch to suppress addition of vertical resistance between bottom of model and groundwater level
      integer   swcfbs             ! Switch for use of coefficient CFBS to convert potential ET into potential E: 0 = no; 1 = yes
      integer   swdiscrvert        ! Switch to convert vertical discretization for soil water quality models: 0 = no; 1 = yes
      ! [GR-BH Task 37] swdivd retired — moved to state%drainage%swdivd
      ! [GR-BH Task 37] swdivdinf retired — moved to state%drainage%swdivdinf
      integer   swdislay           ! Switch to distribute drainage flux vertically with a given position of the top of the model discharge layers: 0 = no; 1 = yes
      integer   swtopdislay(madr)  ! Switch, for each drainage level, to distribute drainage flux vertically with a given position of the top of the model discharge layers: 0 = no; 1 = yes
      integer   swdra              ! Switch for simulation of lateral drainage: 0 = no drainage; 1 = use basic drainage routine; 
                                   !                                            2 = simulate drainage and surface water
      integer   swfrost            ! Switch for reduction of hydraulic conductivity in case of frost: 0 = no; 1 = yes
      integer   swhyst             ! Switch for hysteresis of soil moisture retention function: 0 = no; 1 = yes
      integer   swinco             ! Switch for initial soil moisture condition: 1 = pressure heads; 2 = hydrostatic equilibrium; 
                                   !                                             3 = final pressure heads from previous simulation
      integer   swkmean            ! Switch for mean of hydraulic conductivity: 1 = unweighted arithmic mean, 2 = weighted arithmic mean
                                   !                                            3 = unweighted geometric mean,4 = weighted geometric mean
                                   !                                            5 = unweighted harmonic mean, 6 = weighted harmonic mean
      integer   swkimpl            ! Switch for implicit solution with hydraulic conductivity: 0 = explicit, 1 = implicit
      integer   swliminf           ! Switch for limit of infiltration head to the waterdepth in the channel: 0 = nolimit, 1 = limitation
      integer   swoutputmodflow    ! Switch for extra output file with data for Modflow
      integer   swoxygen           ! Switch oxygen stress: 1 = concept Feddes et al. (1978); 2 = concept Bartholomeus et al. (2008)
      integer   swoxygentype       ! Switch for method oxygen stress calculation: 1 = physical processes; 2 = repro functions
      integer   swpondmx           ! Switch for time dependent maximum amount of ponding (L) on soil surface before runoff starts
      integer   swqhbot            ! Switch for flux-groundwater level relationship: 1 = exponential function; 2 = tabular function
      integer   swcofqhc           ! Switch for additional flux added to exponential flux-groundwater level relationship: 0 = no, 1 = yes
      integer   swredu             ! Switch for reduction of soil evaporation: 0 = no empirical function; 1 = use function of Black; 
                                   !                                           2 = use function of Boesten/Stroosnijder
      integer   swsophy            ! Switch for input of soil hydraulica properties as function parameters (0) or as table (1)
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      integer   ientrytablay(maho,0:matabentries) ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer

      integer   swtopsub           ! Switch for topsoil or subsoil: 1 = topsoil, 2 = subsoil
      
      integer   swrum              ! Switch for RUME output file: 0 = no; 1 = yes
      integer   swini              ! Switch for output files SoilPhysParam.csv and HeatParam.csv: 0 = no; 1 = yes
      integer   swend              ! Switch for output file *.END with end conditions: 0 = no; 1 = end of simulation; 2 = each day
      integer   swwba              ! Switch for output file *.WBA with daily water balance: 0 = no; 1 = yes
      integer   swvap              ! Switch for output file *.VAP with soil profile data (water content, pressure head, 
                                   !        solute concentration, temperature): 0 = no; 1 = yes
      integer   vap                ! Internal number of output file *.VAP with soil profile data 
      integer   wba                ! Internal number of output file *.WBA with cumulative water balance data
      real(8)   aqamp              ! Amplitude of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(8)   aqave              ! Average hydraulic head in deep aquifer (L)
      real(8)   aqper              ! Period of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(8)   aqtmax             ! Time with maximum hydraulic head in deep aquifer (T)
      real(8)   atop(18,6)         ! Parameters of reprofunctions for oxygen stress
      real(8)   basegw             ! Depth of impervious layer (L) for drainage according to Hooghoudt or Ernst
      real(8)   bdens(maho)        ! Array with dry bulk density for each soil layer (M/L3)
      real(8)   c_top(macp)        ! Oxygen concentration at top of compartment(kg/m3)
      
      ! Oxygen stress per-node arrays (scaffolding from OxygenStress subroutine)
      real(8)   o2_d_soil_term1(macp)     ! Pre-calculated soil diffusion term1 per node
      real(8)   o2_d_soil_term2(macp)     ! Pre-calculated soil diffusion term2 per node
      real(8)   o2_gfp100(macp)           ! Gas-filled porosity * 100 per node
      real(8)   o2_capac_term(macp)       ! Water capacity term per node
      real(8)   o2_nmin1(macp)            ! N-1 per node for VG equation
      real(8)   o2_mplus1(macp)           ! M+1 per node for VG equation
      logical   o2_ini_stress             ! O2 stress initialization flag (initialized to .true. via data statement)
      
      real(8)   cfbs               ! Coefficient (-) to convert potential evapotranspiration into potential evaporation
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   cgird              ! Cumulative amount of gross irrigation (L)
      ! real(8)   cinund             ! Cumulative amount of inundation (L)
      ! real(8)   cnird              ! Cumulative amount of net irrigation (L)
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%cofani (GR-BH arc)
      ! real(8)   cofani(maho)       ! Anisotropy coefficient (horizontal / vertical saturated hydraulic conductivity) (-)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%cofgen (ADR 0038)
      ! real(8)   cofgen(21,macp)    ! Array (adjusted by hysteresis) soil hydraulic parameters according to Mualem - van Genuchten for each soil layer
      real(8)   cofqha             ! Coefficient A in exponential relationship between drainage flux and groundwater level (L/T)
      real(8)   cofqhb             ! Coefficient B in exponential relationship between drainage flux and groundwater level (/T)
      real(8)   cofqhc             ! Coefficient C (flux) in exponential relationship between drainage flux and groundwater level (L/T)
      real(8)   cofred             ! Soil evaporation coefficient of Black or Boesten/Stroosnijder
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   cqbot              ! Cumulative amount of water flow through bottom of simulated soil column (L)
      ! real(8)   cqbotdo            ! Cumulative amount of water (L) passed through the soil column bottom in downward direction
      ! real(8)   cqbotup            ! Cumulative amount of water (L) passed through the soil column bottom in upward direction
      ! SS-SWST Phase 2 Task 11 C2: cqdra/cqdrain/cqdrainin/cqdrainout removed — state%surfacewater owns them.
      ! real(8)   cqdra              ! Moved to surfacewater_state_t%cqdra
      ! real(8)   cqdrain(Madr)      ! Moved to surfacewater_state_t%cqdrain(:)
      ! real(8)   cqdrainin(Madr)    ! Moved to surfacewater_state_t%cqdrainin(:)
      ! real(8)   cqdrainout(Madr)   ! Moved to surfacewater_state_t%cqdrainout(:)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   cqprai             ! Cumulative amount of net rain (L)
      ! real(8)   cqssdi             ! Cumulative amount of subsurface drip irrigation (L)
      ! real(8)   cqrot              ! Cumulative amount of extracted water by roots (L)
      ! real(8)   cqtdo              ! Cumulative amount of water (L) passed through the soil surface in downward direction
      ! real(8)   cqtup              ! Cumulative amount of water (L) passed through the soil surface in upward direction
      real(8)   CritDevMasBal      ! Maximum error in water balance (L)
      real(8)   CriterHr           ! Maximum difference of Hroot between iterations; convergence criterium  (L)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   crunoff            ! Cumulative runoff (L)
      ! real(8)   crunon             ! Cumulative amount of runon (L)
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   deepgw             ! hydraulic head in aquifer (L)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%dimoca (ADR 0038)
      ! real(8)   dimoca(macp)       ! Differential soil moisture capacity (/L)
      ! [GR-BH Task 35] disnod(macp+1) retired — moved to state%mesh%disnod
      ! real(8)   drainl(Madr)       ! Moved to drainage_state_t%drainl (ADR 0031)
      real(8)   drares(Madr)       ! Array with drainage resistance (T) for each drainage level
      ! [GR-BH Task 35] dz(macp) retired — moved to state%mesh%dz
      real(8)   dznew(macp)        ! Desired thickness of compartments for soil water quality models (L)
      real(8)   entres             ! Drain entry resistance (T)
      real(8)   es0                !  potential evaporation rate from a wet bare soil [mm/d]
      real(8)   et0                !  potential transpiration rate from a dry crop [mm/d]
      real(8)   ew0                !  potential transpiration rate from a wet crop [mm/d]
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%evp (ADR 0038)
      ! real(8)   evp(macp)          ! Internal evaporation flux of top soil compartments (L/T)
      ! [GR-BH Task 37] FacDpthInf retired — moved to state%drainage%FacDpthInf
      real(8)   ftopdislay(madr)   ! Array with factor for function to determine depth of top of model discharge layer for each drain level, see also swtopdislay (L)
      real(8)   geofac             ! Geometry factor (-) for analytical drainage formula of Ernst
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%gwl (ADR 0038)
      ! real(8)   gwl                ! Groundwater level (L)
      real(8)   gwlconv            ! Maximum difference of groundwater levels between iterations to solve Richards equation
      real(8)   gwli               ! Groundwater level (L) at start of simulation
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   gwlinp             ! Prescribed groundwater level (L) for current time
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%gwlm1 (ADR 0038)
      ! real(8)   gwlm1              ! Groundwater level (L) at former time level
      real(8)   gwltab(mabbc*2)    ! Array with prescribed groundwater level (L) as function of time (T)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%h (ADR 0038)
      ! real(8)   h(macp)            ! Soil water pressure head (L)
      real(8)   h_enpr(macp)       ! Soil water Entry Pressure head for Modified MualemVanGenuchten curve (L)
      real(8)   haqtab(mabbc*2)    ! Array with specified hydraulic head in deep aquifer (L) as function of time (T)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%hatm (ADR 0038)
      ! real(8)   hatm               ! Pressure head of air (L) near the soil surface
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   hbot               ! Soil water pressure head (L) at bottom of soil column
      real(8)   hbotab(mabbc*2)    ! Array with specified pressure head of lowest compartment (L) as function of time (T)
      real(8)   hcomp(macp)        ! Array with prescribed height of numerical compartments (L) for each sublayer
      real(8)   hdrain             ! Mean drainage level (L) to derive regional average groundwater level for bottom boundary condition
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%hleaf (ADR 0036)
      ! real(8)   hleaf              ! Pressure head inside leaves (cm)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%hm1 (ADR 0038)
      ! real(8)   hm1(macp)          ! Soil water pressure head (L) at former time level
      real(8)   hplate             ! Pressure head of ceramic plate below lysimeter
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%hroot (ADR 0036)
      ! real(8)   hroot(macp)        ! Pressure head of a compartment at the root-soil interface (L)
      real(8)   hsublay(macp)      ! Array with prescribed height of sublayers (L)
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   hsurf              ! Soil water pressure head at the soil surface (cm)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%Hxylem (ADR 0036)
      ! real(8)   Hxylem             ! Pressure head in root xylem (L)
      ! real(8) :: igrai         ! Intermediate amount of gross rainfall (L) — [SS-ATM] retired 2026-05-11
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   ies0               ! Intermediate potential evaporation rate from a wet bare soil [cm/d]
      ! real(8)   iet0               ! Intermediate potential transpiration rate from a dry crop [cm/d]
      ! real(8)   iew0               ! Intermediate potential transpiration rate from a wet crop [cm/d]
      ! real(8)   iintc              ! Intermediate amount of rainfall interception by vegetation (L)
      real(8)   infres(Madr)       ! Array with infiltration resistance (T) for each drainage level
      real(8)   inpola(macp)       ! Weight for interpolation between current node and upper node
      real(8)   inpolb(macp)       ! Weight for interpolation between current node and lower node
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%inq (ADR 0038)
      ! real(8)   inq(macp+1)        ! Array with intermediate amounts of water flow between current and upper compartment (L)
      ! SS-SWST Phase 2 Task 11 C2: inqdra/inqdra_in/inqdra_out removed — state%surfacewater owns them.
      ! real(8)   inqdra(Madr,macp)      ! Moved to surfacewater_state_t%inqdra
      ! real(8)   inqdra_in(Madr,macp)   ! Moved to surfacewater_state_t%inqdra_in
      ! real(8)   inqdra_out(Madr,macp)  ! Moved to surfacewater_state_t%inqdra_out
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   inqrot(macp)       ! Array with intermediate amounts of extracted water by roots for each compartment (L)
      ! real(8)   inqssdi(macp)      ! Array with intermediate amounts of subsurface drip irrigation for each compartment (L)
      ! real(8)   qpotrot_day(macp)  ! Array with amounts of potential extracted water by roots for each compartment since start of day (L)
      ! real(8)   qredtot_day(macp)  ! Array with amounts of reduction of extracted water by roots for each compartment since start of day (L)
      ! real(8)   ipondbeg           ! Ponding water layer (L) on soil surface at start of current intermediate period
      ! real(8)   iprec              ! Intermediate amount of gross precipitation + gross irrigation (L)
      ! real(8)   iqbot              ! Intermediate amount of water flow through bottom of simulated soil column (L)
      ! [SS-SWC] iqtdo/iqtup retired — iqinfmax/qinfmax stay legacy (out of scope)
      real(8)   iqinfmax, qinfmax
      ! real(8)   iqdo(macp+1), iqup(macp+1)
      ! SS-SWST Phase 2 Task 11 C2: iqdra removed — state%surfacewater%iqdra owns it.
      ! real(8)   iqdra              ! Moved to surfacewater_state_t%iqdra
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   iqrot              ! Intermediate amount of extracted water by roots (L)
      ! real(8)   iqssdi             ! Intermediate amount of water input via subsurface drip irrigation (L)
      ! real(8)   iqredwet           ! Intermediate amount of reduced root water extraction due to wet conditions (L)
      ! real(8)   iqreddry           ! Intermediate amount of reduced root water extraction due to dry conditions (L)
      ! real(8)   iqredsol           ! Intermediate amount of reduced root water extraction due to salt conditions (L)
      ! real(8)   iqredfrs           ! Intermediate amount of reduced root water extraction due to frost conditions (L)
      ! real(8)   iruno              ! Intermediate amount of runoff (L)
      ! real(8)   irunon             ! Intermediate amount of runon (L)
      real(8)   issnowbeg          ! Amount of snow in soil water equivalent (L) at start of current intermediate period
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%IThetaBeg (ADR 0038)
      ! real(8)   ithetabeg(macp)    ! Array with volumetric soil water contents (-) for each compartment at start of intermediate period
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%k (ADR 0038)
      ! real(8)   k(macp+1)          ! Array with soil hydraulic conductivity (L/T) for each numerical compartment
      real(8)   khbot              ! Horizontal hydraulic conductivity of bottom layer (L/T)
      real(8)   khtop              ! Horizontal hydraulic conductivity of top layer (L/T)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%kmean (ADR 0038)
      ! real(8)   kmean(macp+1)      ! Array with mean soil hydraulic conductivity (L/T) at the interface of current and upper compartment
      real(8)   Kroot              ! Hydraulic, radial conductivity of root tissue (L/T)
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%ksatfit/ksatexm (GR-BH arc)
      ! real(8)   ksatfit(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: fitted on VG based on lab data
      ! real(8)   ksatexm(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: examined in lab or field
      real(8)   ksatthr(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: to interpolate VG and Ksatexm
      real(8)   kstem              ! Conductance in the path from leaf to root xylem (/d)
      real(8)   kvbot              ! Vertical hydraulic conductivity of bottom layer (L/T)
      real(8)   kvtop              ! Vertical hydraulic conductivity of top layer (L/T)
      ! [GR-BH Task 37] L(Madr) retired — moved to state%drainage%L
      ! real(8) :: ldwet         ! Length of dry period (L) as used in Black's model — [SS-ATM] retired 2026-05-11
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%mfluxtable (ADR 0036)
      ! real(8)   mfluxtable(maho,801)  ! Reference table with matric flux potential of each soil layer (L2/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%mflux (ADR 0036)
      ! real(8)   mflux(macp)        ! Actual matric flux potential of each node (L2/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%mroot (ADR 0036)
      ! real(8)   mroot(macp)        ! Matrix flux head of a compartment at the root-soil interface (L2/T)
      real(8)   OxygenIntercept(6) ! Parameters of reproduction function for oxygen stress according to Bartholomeus
      real(8)   OxygenSlope(6)     ! Parameters of reproduction function for oxygen stress according to Bartholomeus
      real(8)   paramvg(21,maho)   ! Array with input values of soil hydraulic parameters according to Mualem - van Genuchten for each soil layer
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   pegwl              ! Perched groundwater level (L)
      ! real(8)   pond               ! Height of ponding layer (L)
      ! real(8)   pondini            ! Ponding water layer (L) on soil surface at start of current water balance period
      ! real(8)   pondm1             ! Ponding water layer (L) on soil surface at former time level
      real(8)   pondmx             ! Maximum amount of ponding (L) on soil surface before runoff starts
      real(8)   pondmxtab(2*mairg) ! Table with time-dependent input (date,value) for maximum amount of ponding (L) on soil surface before runoff starts
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%q (ADR 0038)
      ! real(8)   q(macp+1)          ! Soil water flux between current compartment and upper compartment (L/T)
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   qbot               ! Water flux through bottom of simulated soil column (L/T)
      real(8)   qbotab(mabbc*2)    ! Array with specified bottom flux (L/T) as function of time (T)
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   qbot_nonfrozen     ! Water flux through bottom of non-frozen soil column (L/T)
      ! real(8)   qdra(Madr,macp)    ! Moved to drainage_state_t%qdra (ADR 0031)
      ! real(8)   qdrain(Madr)       ! Moved to drainage_state_t%qdrain (ADR 0031)
      real(8)   qdraincomp(macp)   ! Total lateral drainage flux (L/T) for each compartment
      real(8)   qdrtab(50)         ! Array with lateral drainage flux (L/T) as function of groundwater level (L)
      ! SS-SWST Phase 2 Task 11 C2: qdrtot removed — state%surfacewater%qdrtot owns it.
      ! real(8)   qdrtot             ! Moved to surfacewater_state_t%qdrtot
      real(8)   qimmob(macp)       ! Soil water flux between mobile and immobile fraction in case of fingered flow (L/T)
      real(8)   qssdisum           ! Total subsurface irrigation flux (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qrosum (ADR 0036)
      ! real(8)   qrosum             ! Total root water extraction flux (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qredwetsum (ADR 0036)
      ! real(8)   qredwetsum         ! Total reduction of root water extraction due to wet conditions (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qreddrysum (ADR 0036)
      ! real(8)   qreddrysum         ! Total reduction of root water extraction due to dry conditions (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qredsolsum (ADR 0036)
      ! real(8)   qredsolsum         ! Total reduction of root water extraction due to salt conditions (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qredfrssum (ADR 0036)
      ! real(8)   qredfrssum         ! Total reduction of root water extraction due to frost conditions (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qrot (ADR 0036)
      ! real(8)   qrot(macp)         ! Array with root water extraction flux for each compartment (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qpotrot (ADR 0036)
      ! real(8)   qpotrot(macp)      ! Array with potential root water extraction flux for each compartment (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qredwet (ADR 0036)
      ! real(8)   qredwet(macp)      ! Array with reduction of root water extraction due to wet conditions for each compartment (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qreddry (ADR 0036)
      ! real(8)   qreddry(macp)      ! Array with reduction of root water extraction due to dry conditions for each compartment (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qredsol (ADR 0036)
      ! real(8)   qredsol(macp)      ! Array with reduction of root water extraction due to salt conditions for each compartment (L/T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qredfrs (ADR 0036)
      ! real(8)   qredfrs(macp)      ! Array with reduction of root water extraction due to frost conditions for each compartment (L/T)
      real(8)   qssdi(macp)        ! Array with water input via subsurface drip irrigation for each compartment (L/T)
      real(8)   dt_SSDI_event      ! Length of SSDI irrigation event (T)
      
      ! SSDI persistent state (moved from irrigation.f90 local SAVE)
      integer   swssdi_irr         ! Switch: SSDI active (0=no, 1=yes)
      integer   nod_ssdi_irr(2)    ! Upper and lower nodes for SSDI
      integer   ssdi_schedule_irr  ! Schedule type (0=fixed dates, 1=internal)
      integer   ssdi_sched_type_irr ! Internal schedule type (1=Tact/Tpot, 2=h, 3=theta)
      integer   nod_ssdi_sensor_irr ! Sensor node (if ssdi_sched_type > 1)
      real(8)   ssdi_threshold_irr ! Threshold value for scheduling
      real(8)   ssdi_threshold_z_irr ! Depth for threshold value
      real(8)   ssdi_amount_irr    ! Amount of scheduled irrigation (cm)
      real(8)   ssdi_appl_rate_irr ! Application rate (cm/d)
      integer   sw_interval_irr    ! Switch for minimum interval
      integer   days_interval_irr  ! Minimum days between applications
      integer   days_counter_irr   ! Days since previous application
      integer   nirri_ssdi_irr     ! SSDI counter/entry point
      real(8)   ssdi_date_irr(mairg)   ! Fixed irrigation dates
      real(8)   ssdi_rate_f_irr(mairg) ! Fixed irrigation rates (cm/d)
      real(8)   ssdi_amount_f_irr(mairg) ! Fixed irrigation amounts (cm)
      
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   qtop               ! Water flux through soil surface (L/T)
      real(8)   relsatthr(maho)    ! Array with relative saturation (-) for each soil layer: to interpolate VG and Ksatexm
      real(8)   ResultsOxygenStress(19,macp) ! array with results for OxygenStress; for output only 
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   reva               ! Actual soil evaporation rate (L/T)
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%rfcp (ADR 0034)
      ! real(8)   rfcp(macp)         ! Reduction factor for frozen conditions in each model compartment (-)
      real(8)   rimlay             ! Vertical resistance of aquitard (T)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%rmax (ADR 0036)
      ! real(8)   rmax(macp)         ! Radius around roots in which water is extracted (L)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%rootphi (ADR 0036)
      ! real(8)   RootPhi(macp)      ! Factor Phi of a compartment used in drought reduction of De Jong van Lier et al. (T/L)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%rootrho (ADR 0036)
      ! real(8)   RootRho(macp)      ! Factor Rho of a compartment used in drought reduction of De Jong van Lier et al. (/L2)
      real(8)   rsigni             ! Minimum amount of rainfall (L) which resets the empirical soil evaporation reduction models
      real(8)   rsoil              ! Soil resistance of wet soil of PMdirect (T/L)
      real(8)   rsro               ! Drainage resistance for surface runoff (T)
      real(8)   rsroexp            ! Exponent to calculate surface runoff (T)
      integer   swuseCN            ! Switch for usage of Curve Number method for runoff [0 (default),1]
      integer   wc_cor
      real(8)   CNref              ! Input value for Curve Number (CN) [0.01,100]
      real(8)   CNrefTAB(mayrs*5), CNtimTAB(mayrs*5)
      integer   iCNtab
      ! [SS-SWC] retired 2026-05-12 — crunoffCN/irunoCN moved to state%soilwater (ADR 0038)
      real(8)   CNdry, CNwet, ThetaRef, wc10, Runoff_CN  ! crunoffCN, irunoCN now in state%soilwater
      real(8)   Rxylem             ! Mean radius of xylem tube inside roots (L)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%runon (ADR 0038)
      ! real(8)   runon              ! Water runon flux (L/T)
      real(8)   runonarr(maday)    ! Array with runon (L) data for each day
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   runots             ! Amount of runoff during a time step (L)
      ! real(8) :: saev          ! Cumulative actual evaporation (L) Boesten/Stroosnijder — [SS-ATM] retired 2026-05-11
      real(8)   shape              ! Shape factor: ratio between the mean and the maximum groundwater level elevation above the drainage base (-)
      real(8)   sinamp             ! Amplitude of prescribed bottom flux (L/T) in case of sine function
      real(8)   sinave             ! Average value of prescribed bottom flux (L/T) in case of sine function
      real(8)   sinmax             ! Time of the year with maximum bottom flux in case of prescribed sine function
      ! real(8) :: spev          ! Cumulative potential evaporation (L) Boesten/Stroosnijder — [SS-ATM] retired 2026-05-11
      real(8)   sptab(7,macp,matab)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      real(8)   sptablay(7,maho,matab) ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer
      real(8)   StepHr             ! Maximum difference of Hroot and Hxylem between iterations; convergence criterium  (L)
      real(8)   taccur             ! Maximum absolute difference between simulated and calculated potential transpiration rate (cm/d)
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%Tactual (ADR 0036)
      ! real(8)   Tactual            ! Actual transpiration at former iteration in JongvanLier (cm/d)
      real(8)   tau                ! Minimum pressure head difference (L) to change from wetting to drying in case of hysteresis
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   theta(macp)        ! Volumic soil water content (-)
      ! real(8)   thetar(macp)       ! Residual volumic soil water content (-) for each numerical compartment
      ! real(8)   thetas(macp)       ! Saturated volumic soil water content (-) for each numerical compartment
      ! real(8)   thetm1(macp)       ! Volumic soil water content (-) at former time level
      ! real(8)   thetsl(maho)       ! Saturated volumic water content (-) for each soil layer
      real(8)   twilt(macp)        ! Pressure head of a compartment at wilting point (L)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   volact             ! Water storage (L) of soil column at current time level
      ! real(8)   volini             ! Water storage (L) of soil column at start of simulation
      ! real(8)   volm1              ! Water storage (L) of soil column at former time level
      ! SS-SWST Phase 2 Task 11 C2: vtair removed — state%surfacewater%vtair owns it.
      ! real(8)   vtair              ! Moved to surfacewater_state_t%vtair
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%wbalance (ADR 0038)
      ! real(8)   wbalance           ! Cumulative water balance error (L)
      ! real(8)   wetper(Madr)       ! Moved to drainage_state_t%wetper (ADR 0031)
      ! [GR-BH Task 35] z(macp) retired — moved to state%mesh%z
      ! [GR-BH Task 35] ztopcp(macp) retired — moved to state%mesh%ztopcp
      ! [GR-BH Task 35] zbotcp(macp) retired — moved to state%mesh%zbotcp
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%zfrostbot (ADR 0034)
      ! real(8)   zfrostbot          ! Depth of bottom of frost layer (L)
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%zfrosttop (ADR 0034)
      ! real(8)   zfrosttop          ! Depth of top of frost layer (L)
      ! [GR-BH Task 37] zbotdr(Madr) retired — moved to state%drainage%zbotdr
      real(8)   zi(macp)           ! Array with soil depths (L) used to specify initial soil water pressure heads
      real(8)   zintf              ! Depth (L) at which fine top layer ends and coarse sub layer starts
      ! real(8)   ztopdislay(Madr)   ! Moved to drainage_state_t%ztopdislay (ADR 0031)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flDrain (ADR 0041)
      ! logical   fldrain            ! moved to state%timecontrol%flDrain
      logical   FlHydrLift         ! Flag indicating release of water from root to soil is allowed
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%fllowgwl (ADR 0038)
      ! logical   fllowgwl           ! Flag indicating precribed groundwater level below bottom soil column
      logical   flrunon            ! Flag indicating the existance of runon
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! logical   ftoph              ! Flag indicating that the pressure head is prescribed at the soil surface
      character(len=16) drfil      ! Name of drainage input file
      character(len=80) pathdrain  ! Path to folder with drainage input files

! --- heat variables
      integer   nheat              ! Number of initial soil temperatures as provided in the input
      integer   swbotbhea          ! Switch for bottom boundary condition: 1 = heat flux is zero; 2 = prescribed temperature
      integer   swtopbhea          ! Switch for top boundary condition: 1 = use air temperatures; 2 = read measured surface temperatures
      integer   swcalt             ! Switch for method of soil water heat flow simulation: 1 = analytical method; 2 = numerical method
      integer   swhea              ! Switch for simulation of soil heat flow: 0 = no; 1 = yes
      integer   swtem              ! Switch for output file *.TEM with soil temperatures: 0 = no; 1 = yes
      integer   tem                ! Internal number of output file *.TEM with soil temperatures
      real(8)   ddamp              ! Damping depth (L) of temperature wave in soil
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%fclay (ADR 0034)
      ! real(8)   fclay(macp)        ! Array with gravimetric content of clay (g/g mineral parts) of each numerical compartment
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%forg (ADR 0034)
      ! real(8)   forg(macp)         ! Array with gravimetric organic matter content (g/g mineral parts) of each numerical compartment
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%fquartz (ADR 0034)
      ! real(8)   fquartz(macp)      ! Array with gravimetric content of sand+silt (g/g mineral parts) of each numerical compartment
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%orgmat/pclay/psand/psilt (GR-BH arc)
      ! real(8)   orgmat(maho)       ! Array with gravimetric organic matter content (g/g mineral parts) for each soil layer
      ! real(8)   pclay(maho)        ! Array with gravimetric clay content (g/g mineral parts) for each soil layer
      ! real(8)   psand(maho)        ! Array with gravimetric sand content (g/g mineral parts) for each soil layer
      ! real(8)   psilt(maho)        ! Array with gravimetric silt content (g/g mineral parts) for each soil layer
      real(8)   tampli             ! Amplitude of prescribed annual temperature wave (�C) at soil surface
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%tebot (ADR 0034)
      ! real(8)   tebot              ! Temperatures (�C) at bottom of soil profile
      real(8)   tembtab(mabbc*2)   ! Array with specified bottom temperature (�C) as function of time (T)
      real(8)   temtoptab(mabbc*2) ! Array with specified soil surface temperature (�C) as function of time (T)
      real(8)   tfroststa          ! Soil temperature (�C) where reduction of water fluxes starts
      real(8)   tfrostend          ! Soil temperature (�C) where reduction of water fluxes ends
      real(8)   timref             ! Time in the year (T) with top of prescribed sine temperature wave
      real(8)   tmean              ! Prescribed mean annual temperature (�C) at soil surface
      ! [SS-HEAT] NOTE: tsoil retained as config-staging buffer for temperature.f90 task=1 init.
      ! The runtime soil-temperature compute state is state%heat%tsoil(numnod). (ADR 0034)
      real(8)   tsoil(macp)        ! Config-staging: nheat initial temperature profile entries (�C) for each compartment
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%tetop (ADR 0034)
      ! real(8)   tetop              ! Temperatures (�C) at top of soil profile (under snow cover)
      real(8)   zh(macp)           ! Array with soil depths (L) used to specify initial soil temperatures
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%heacap (ADR 0034)
      ! real(8)   heacap(macp)       ! Array with heat capacity for all compartments (J/cm3/K)
      ! [SS-HEAT] retired 2026-05-10 — moved to state%heat%heacon (ADR 0034)
      ! real(8)   heacon(macp)       ! Array with heat conductivity for all compartments (J/cm/K/d)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flTemperature (ADR 0041)
      ! logical   fltemperature      ! moved to state%timecontrol%flTemperature

! --- snow variables
      integer   snw                ! Internal number of output file *.SNW with snow pack data
      integer   swsnow             ! Switch for simulation of snow accumulation and melt: 0 = no; 1 = yes
      integer   swsublim           ! Switch for suppressing simulation of sublimation of snow: 1 = suppress ! Adaptation 3 for PEARL-MACRO
      ! ========================================================================
      ! [SS-ATM] retired 2026-05-11 — snow scalars migrated to state%atmosphere (flat) / %cumu / %intr
      !   See: ADR 0037 (A-2.7)
      ! ========================================================================
      ! real(8) :: cgsnow    ! Cumulative amount of gross snow fall (L water)
      ! real(8) :: cmelt     ! Cumulative amount of melted snow (L water)
      ! real(8) :: csnrai    ! Cumulative amount of net snow fall (L water)
      ! real(8) :: csubl     ! Cumulative amount of sublimated snow (L water)
      ! real(8) :: gsnow     ! Gross snow rate (L/T)
      ! real(8) :: igsnow    ! Incremental amount of gross snow fall (L water)
      ! real(8) :: isnrai    ! Incremental amount of net snow fall (L water)
      ! real(8) :: isubl     ! Incremental amount of sublimated snow (L water)
      ! real(8) :: melt      ! Melting rate (L/T)
      ! real(8) :: slw       ! liquid water stored in snow pack
      real(8)   snowcoef           ! Snow melt factor (-)
      ! real(8) :: snowinco  ! Amount of snow (L water) at start of balance period
      ! real(8) :: snrai     ! Net rain rate on snow pack (L/T)
      ! real(8) :: ssnow     ! Amount of snow (L water)
      ! real(8) :: subl      ! Sublimation rate (L/T)
      real(8)   TePrRain           ! Temperature above which all precipitation is rain,[ 0.0...5.0 oC, R]
      real(8)   TePrSnow           ! Temperature below which all precipitation is snow,[-5.0...0.0 oC, R]
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flSnow (ADR 0041)
      ! logical   flsnow             ! moved to state%timecontrol%flSnow

! --- solute variables
      integer   nconc              ! Number of initial solute concentrations as provided in the input
      integer   sba                ! Internal number of output file *.SBA with cumulative solute balance components
      integer   swbr               ! Switch to consider mixed reservoir for solute breakthrough in the saturated zone: 0 = no; 1 = yes
      integer   swbotbc            ! Switch for bottom boundary condition of solute-concentration (see *.SWP input file for overview)
      integer   swsolu             ! Switch for simulation of solute transport: 0 = no; 1 = yes
      integer   swsp               ! Switch (in case of solute transport) for simulation of sorption  0 = no; 1 = yes
      real(8)   AgeGwl1m           ! Age (d) of groundwater in upper 1 meter of saturated zone
      real(8)   bexp               ! Exponent in decomposition reduction factor due to dryness (-)
      real(8)   cdrain             ! Mean solute concentration in aquifer or drainage system (M/L3 water)  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(8)   cirr               ! Solute concentration (M/L3) in irrigation water
      real(8)   cml(macp)          ! Array with solute concentration (M/L3 water) in mobile region
      real(8)   cmsy(macp)         ! Array with dissolved + adsorbed solute concentration (M/L3 soil volume) in mobile region
! real(8)   cpond              ! Moved to solute_state_t%cpond (ADR 0032)
      real(8)   cpre               ! Solute concentration (M/L3) in precipitation
      real(8)   cref               ! Reference solute concentration (M/L3) for Freundlich adsorption
! real(8)   cseep              ! Moved to solute_state_t%cseep (ADR 0032)
      real(8)   cseeptab(mabbc*2)  ! Array with Mean solute concentration in upward seepage water at bottom of profile (M/L3 water) as function of time (T)
! real(8)   csurf              ! Moved to solute_state_t%csurf (ADR 0032)
      real(8)   daquif             ! Thickness of saturated aquifer (L) to calculate solute breakthrough to surface water
      real(8)   ddif               ! Molecular diffusion coefficient (L2/T)
      real(8)   decpot(maho)       ! Array with Potential decomposition rate (/T) for each soil layer
      real(8)   decsat             ! Decomposition rate in aquifer (/T)
! real(8)   dectot             ! Moved to solute_state_t%dectot (ADR 0032)
! real(8)   imdectot           ! Moved to solute_state_t%imdectot (ADR 0032)
      real(8)   dtsolu             ! Maximum time step (T) for accurate numerical solution of solute transport equation  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(8)   fdepth(maho)       ! Array with reduction factor for decomposition (-) for each soil layer
      real(8)   frexp              ! Array with Freundlich exponent (-) for solute adsorption
      real(8)   gampar             ! Reduction factor for decomposition due to low temperatures (/C)
      real(8)   icAgeBot           ! Incremental (over output interval) age (d) of groundwater leaving bottom comp.
      real(8)   icAgeDra(madr)     ! Incremental (over output interval) age (d) of groundwater leaving bottom comp.
      real(8)   icAgeRot           ! Incremental (over output interval) age (d) of groundwater leaving by root uptake
      real(8)   icAgeSur           ! Incremental (over output interval) age (d) of groundwater leaving by surface runoff
      real(8)   isqbot             ! Solute flux at the bottom of the soil column (M/L2/T)  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(8)   isqtop             ! Solute flux through the soil top surface (M/L2/T)  [AgeTracer dead-code dep, keep until agetracer_state_t]
      real(8)   kf(maho)           ! Array with Freundlich coefficient (L3/M) for solute adsorption for each soil layer
      real(8)   kfsat              ! Linear adsorption coefficient in aquifer (L3/M)
      real(8)   ldis(maho)         ! Array with Solute dispersion length (L) for each soil layer
      real(8)   poros              ! Porosity of aquifer (-) to calculate solute breakthrough
      real(8)   rottot             ! Cumulative amount of solutes (M/L2) extracted by plant roots  [AgeTracer dead-code dep, keep until agetracer_state_t]
! real(8)   imrottot           ! Moved to solute_state_t%imrottot (ADR 0032)
      real(8)   rtheta             ! Minimum volumetric water content (-) for potential decomposition
      real(8)   salthead           ! Conversion salt concentration (mg/cm3) into osmotic head (cm) [0..1000.0 cm/(mg/cm3), R]
      real(8)   saltmax            ! Threshold salt concentration in soil water  [0..100 mg/cm3, R]
      real(8)   saltslope          ! Decline of rootwater uptake above threshold [0..1.0 cm3/mg, R]
! real(8)   samcra             ! Moved to solute_state_t%samcra (ADR 0032)
      real(8)   samini             ! Total amount of solutes (M/L2) in soil profile at start of current balance period  [AgeTracer dead-code dep, keep until agetracer_state_t]
! real(8)   sampro             ! Moved to solute_state_t%sampro (ADR 0032)
! real(8)   solbal             ! Moved to solute_state_t%solbal (ADR 0032)
! real(8)   sqbot              ! Moved to solute_state_t%sqbot (ADR 0032)
! real(8)   imsqbot            ! Moved to solute_state_t%imsqbot (ADR 0032)
      real(8)   sqdra              ! Total amount of solutes (M/L2) transported to drainage canals  [AgeTracer dead-code dep, keep until agetracer_state_t]
! real(8)   imsqdra            ! Moved to solute_state_t%imsqdra (ADR 0032)
! real(8)   sqirrig            ! Moved to solute_state_t%sqirrig (ADR 0032)
! real(8)   imsqirrig          ! Moved to solute_state_t%imsqirrig (ADR 0032)
! real(8)   sqprec             ! Moved to solute_state_t%sqprec (ADR 0032)
! real(8)   imsqprec           ! Moved to solute_state_t%imsqprec (ADR 0032)
! real(8)   sqrap              ! Moved to solute_state_t%sqrap (ADR 0032)
! real(8)   sqsur              ! Moved to solute_state_t%sqsur (ADR 0032)
      real(8)   tscf               ! Relative uptake of solutes by roots (-)
      real(8)   zc(macp)           ! Array with soil depths (L) used to specify initial solute concentrations
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flSolute (ADR 0041)
      ! logical   flsolute           ! moved to state%timecontrol%flSolute
      logical   flAgeTracer        ! Flag indicating simulation of Ageing (groundwater age)

! --- age tracer state variables (moved from local SAVE in AgeTracer subroutine)
      real(8)   Ageirr             ! Age of irrigation water (d)
      real(8)   Agedrain           ! Age of drainage water (d)
      real(8)   Agepre             ! Age of precipitation (d)
      real(8)   Agepond            ! Age of ponding water (d)
      real(8)   Agepondm1          ! Age of ponding water previous timestep (d)
      real(8)   icAgetopupw        ! Incremental age leaving top compartment upward (d)
      real(8)   icAgetopdwn        ! Incremental age entering top compartment downward (d)
      ! [MACRO-RETIRE 2026-05-12] ADR 0040 — macropore subsystem retired from
      ! the rescue branch. The following globals are kept as permanent-zero
      ! placeholders because cross-subsystem code reads them under dead
      ! `if (flMacroPore)` branches (flMacroPore is forced .false. in init).
      ! See legacy/swap-4.2.0 for the original SWAP 4.2.0 implementation.
      ! [GR-BH Task 36] retired 2026-05-13 — ArMpSs made local in soilhydraulics/solute/agetracer (GR-BH arc, ADR 0040 complete)
      ! real(8)   ArMpSs             ! Area fraction of macropores at soil surface (-)  [retired-zero]

! --- macropore Input parameters (all retired-zero — see [MACRO-RETIRE 2026-05-12])
      integer SwSoilShr(MaHo)      ! [retired-zero] kept: soilhydraulics shrinkage branch
      real(8) ThetCrMp(MaHo)       ! [retired-zero] kept: soilhydraulics shrinkage branch
      real(8) Z_Tp                 ! [retired-zero] kept: ArMpTp/ArMpSs gating
      real(8) CritUndSatVol        ! [retired-zero] kept: waterbalance watertable() arg
! --- macropore variables (selected — most retired by deletion)
      real(8) ArMpTp               ! [retired-zero] kept: ArMpSs assignment
      real(8) cQMpLatSs            ! [retired-zero] kept: soilhydraulics zero-write
      real(8) cQMpOutDrRap         ! [retired-zero] kept: waterbalance wbalance term
      real(8) dFdhMp(MaCp)         ! [retired-zero] kept: soilhydraulics dFdhM term (always 0)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%dtold (ADR 0041)
      ! real(8) dtold                ! moved to state%timecontrol%dtold
      real(8) DiPoCp(MaCp)         ! [retired-zero] kept: soilgrid refinement
      real(8) iQMpOutDrRap         ! [retired-zero] kept: swap_csv/swapoutput DRAINAGE accumulator
      real(8) iQInTopLatDm1        ! [retired-zero] kept: waterbalance .BMA writer
      real(8) iQInTopLatDm2        ! [retired-zero] kept: waterbalance .BMA writer
      real(8) iQInTopVrtDm1        ! [retired-zero] kept: waterbalance .BMA writer
      real(8) iQInTopVrtDm2        ! [retired-zero] kept: waterbalance .BMA writer
      real(8) IAvFrMpWlWtDm1(MaCp) ! [retired-zero] kept: soilgrid refinement
      real(8) IAvFrMpWlWtDm2(MaCp) ! [retired-zero] kept: soilgrid refinement
      real(8) iQExcMtxDm1Cp(MaCp)  ! [retired-zero] kept: soilgrid refinement
      real(8) iQExcMtxDm2Cp(MaCp)  ! [retired-zero] kept: soilgrid refinement
      real(8) iQOutDrRapCp(MaCp)   ! [retired-zero] kept: soilgrid refinement
      real(8) IWaSrDm1Beg          ! [retired-zero] kept: waterbalance .BMA writer
      real(8) IWaSrDm2Beg          ! [retired-zero] kept: waterbalance .BMA writer
      real(8) WaSrDm1              ! [retired-zero] kept: waterbalance .BMA writer
      real(8) WaSrDm2              ! [retired-zero] kept: waterbalance .BMA writer
      real(8) WaSrDm1Ini           ! [retired-zero] kept: waterbalance wbalance term
      real(8) WaSrDm2Ini           ! [retired-zero] kept: waterbalance wbalance term
      real(8) VlMpStDm1(MaCp)      ! [retired-zero] kept: soilgrid refinement
      real(8) VlMpStDm2(MaCp)      ! [retired-zero] kept: soilgrid refinement
      integer IcTopMP              ! [retired-zero] kept: cross-subsystem ref
      integer IDecMpRat            ! [retired-zero] kept: soilhydraulics convergence counter
      real(8) QExcMpMtx(MaCp)      ! [retired-zero] kept: waterbalance use clause
      real(8) QMaPo                ! [retired-zero] kept: waterbalance qbot term
      real(8) QRapDra              ! [retired-zero] kept: surfacewater drainage terms
      ! [MACRO-RETIRE 2026-05-12] note: NumLevRapDra/RapDraReaExp/RapDraResRef
      ! belong to the drainage subsystem (set by drainage_config), NOT
      ! macropore — they were grouped here legacy-style. Keep them.
      integer NumLevRapDra         ! Number of drainage levels for rapid drainage (drainage feature)
      real(8) RapDraReaExp         ! Reaction coefficient for rapid drainage (drainage feature)
      real(8) RapDraResRef(Madr)   ! Reference rapid drainage resistance (drainage feature)
      logical FlDecMpRat           ! [retired-zero] kept: soilhydraulics convergence sentinel
      logical flmacropore          ! [retired-zero] forced .false. in init — guards dead branches
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! [MACRO-RETIRE 2026-05-12] All other macropore globals retired.
      ! See ADR 0040 and legacy/swap-4.2.0 for the original SWAP 4.2.0
      ! macropore implementation. The retired set includes ~24 cumulative
      ! accumulators (cQMp*, iQMp*, iQExcMtx*, iQInTop*, IWaSrDm*,
      ! IAvFrMpWlWt*), area/depth parameters (NumSbDm, SwDarcy, SwDrRap,
      ! SwPowM, SwShrInp, SwSorp, CritUndSatVol, DiPoMa, DiPoMi, FKcovlay,
      ! GeomFac, PndmxMp, PowM, PpIcSs, RapDra*, Rzah, ShapeFacMp, Sorp*,
      ! ShrPar[A-E], Spoint, VlMpStSs, Z_Ah, Z_Ic, Z_MB50, Z_St, ZDiPoMa,
      ! ZnCrAr, KsMpSs, KsatCovLay, PpDmCp, PpIcTpMp, ArMpSs and 27 more)
      ! and macropore work arrays (ICpBt*, ICpTp*, ICpSat*, NnCrAr,
      ! ArMpTpDm, AwlCorFac, FrMpWalWet, KDCrRlRef, QExcMtxDmCp,
      ! QIn*SatDmCp, QInTop*Dm, QOut*DmCp, SorpDmCp, ThtSrpRefDmCp,
      ! TimAbsCumDmCp, VlMpDm*, WaSrMp*, WaSr/WaLev/VlMpDm1/2 etc., ZBtDm,
      ! ZWaLevDm, flBegin, flDraTub, FlEndSrpEvt, IcTopMP, IDecMpRat,
      ! NumDm, NumLevRapDra, SwBma, SubsidCp, VlMpDyCp, VlMpStCp,
      ! VlMpStDm1/2, dFdhMp is kept retired-zero, etc.).

! --- surface water variables
      integer swswb,swdrf,swsrf,swallo(Madr),swdtyp(Madr)
      ! [GR-BH Task 37] swnrsrf retired — moved to state%drainage%swnrsrf
      integer swqhr,swsec,nrpri,nrsec,nmper,swman(mamp)
      ! [GR-BH Task 37] SwTopnrsrf retired — moved to state%drainage%swtopnrsrf
      ! SS-SWST Phase 2 Task 11 C2: numadj/imper removed — state%surfacewater owns them.
      integer nqh(mamp),drf,swb,nphase(mamp),nodhd(mamp)
      ! numadj removed (surfacewater_state_t%numadj)
      integer intwl(mamp)
      ! imper removed (surfacewater_state_t%imper)
      integer nowltab(madr)
      real(8) widthr(Madr),taludr(Madr),rdrain(Madr),rsurfdeep
      real(8) rsurfshallow,rinfi(Madr),rentry(Madr),rexit(Madr)
      real(8) gwlinf(Madr),wlptab(2*mawlp)
      real(8) impend(mamp)
      real(8) wldip(mamp),wscap(mamp),hbweir(mamp)
      ! SS-SWST Phase 2 Task 11 C2: wlstar removed — state%surfacewater%wlstar owns it.
      real(8) osswlm,wlp,alphaw(mamp),betaw(mamp)
      ! wlstar removed (surfacewater_state_t%wlstar)
      real(8) wls1_init    ! TOML pipeline: initial wls1 = wlact - altcu; altcu=0 is enforced by drainage_config_validate
      real(8) dropr(mamp*mamte),hdepth(mamp*mamte)
      real(8) gwlcrit(mamp,mamte),hcrit(mamp,mamte),vcrit(mamp,mamte)
      real(8) hqhtab(mamp,mamte),qqhtab(mamp,mamte)
      real(8) wlsman(mamp,mamte)
      ! SS-SWST Phase 2 Task 11 C2: sttab removed — state%surfacewater%sttab owns it.
      real(8) wlstab(2*mawls)
      ! sttab(22,2) removed (surfacewater_state_t%sttab)
      ! SS-SWST Phase 2 Task 11 C2: swstini/swst/wlsbak removed — state%surfacewater owns them.
      ! real(8) swstini   ! Moved to surfacewater_state_t%swstini
      ! real(8) swst      ! Moved to surfacewater_state_t%swst
      ! real(8) wlsbak(4) ! Moved to surfacewater_state_t%wlsbak
      ! SS-SWST Phase 2 Task 11 C2: cqdrd/cwsupp/cwout/wls removed — state%surfacewater owns them.
      real(8) cofintfl,expintfl
      ! cqdrd removed (surfacewater_state_t%cqdrd)
      ! cwsupp removed (surfacewater_state_t%cwsupp)
      ! cwout removed  (surfacewater_state_t%cwout)
      ! wls removed    (surfacewater_state_t%wls)
      ! SS-SWST Phase 2 Task 11 C2: hwlman/wlsold removed — state%surfacewater owns them.
      real(8) owltab(Madr,2*maowl)
      ! real(8) qdrd                  ! Moved to drainage_state_t%qdrd (ADR 0031)
      ! hwlman removed (surfacewater_state_t%hwlman)
      ! wlsold removed (surfacewater_state_t%wlsold)
      ! SS-SWST Phase 2 Task 11 C2: overfl removed — state%surfacewater%overfl owns it.
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flSurfaceWater (ADR 0041)
      ! logical flsurfacewater
      ! overfl removed (surfacewater_state_t%overfl)

      ! Preparation before crop growth
      logical   flCropPrep         ! Flag indicating if ploughing opportunity has been realized
      real(8)   zPrep              ! z-level for monitoring work-ability for the crop     
      real(8)   hPrep              ! maximum pressure head during preparation
      integer   MaxPrepDelay       ! maximum delay of preparation (starting from begin of growing season)
      integer   PrepDelay          ! delay of preparation
      real(8)   dhPrep             ! overshoot of pressure head for work-ability during preparation
      
      ! Sowing before crop growth
      logical   flCropSow          ! Flag indicating if sowing opportunity has been realized
      real(8)   zSow               ! z-level for monitoring work-ability for the crop
      real(8)   hSow               ! maximum pressure head during sowing
      real(8)   zTempSow           ! z-level for monitoring temperature for sowing   
      integer   MaxSowDelay        ! maximum delay of sowing (starting from begin of growing season)
      integer   SowDelay           ! delay of delay
      real(8)   TempSow            ! temperature for sowing   
      real(8)   dhSow              ! overshoot of pressure head for work-ability during sowing
      real(8)   dtempSow           ! undershoot of temperature for sowing at end of available period
      
      ! Germination before crop growth
      logical   flCropGerm         ! Flag indicating if germination has been realized
      real(8)   zgerm              ! z-level for monitoring temperature for germination
      integer   DayGerm            ! Daynumber of germiniation

      ! Harvest of crop growth
      integer   swharv             ! switch for timing of harvest (0=cropend, 1=maturity)
      real(8)   hmow               ! maximum pressure head during mowing
      real(8)   hgrz               ! maximum pressure head during grazing
      real(8)   zmow               ! z-level for monitoring wkability for the crop during mowing
      real(8)   zgrz               ! z-level for monitoring wkability for the crop during start of grazing

! --- root development of dynamic crop growth and oxygen stress
      integer   swWrtNonox         ! switch for checking oxygen stress of root zone development
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%flWrtNonox (ADR 0036)
      ! logical   flWrtNonox         ! Flag indicating whether root development is retatarded by oxygen stress 
      real(8)   aeratecrit         ! threshold to stop root zone development in case of oxygenstress; 0.0 minimum oxygen stress, 1.0 maximum oxygenstress [0.0001..1.0 -, R]

      ! Initialize o2_ini_stress to .true. (needed for first call to OxygenStress)
      data o2_ini_stress /.true./

      end module variables
