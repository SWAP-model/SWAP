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
      ! [GR-ATM 2026-05-23] logf retired — swap_log owns the log-file unit

      ! [SS-GR-FINAL D1] ex_tlast retired — handle_exchange dead; 0 consumers
      
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
      ! [SS-GR-FINAL D1] o2_initialized retired — 0 consumers (o2_ini_stress is the one with consumers)
      
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
      ! [GR-ATM 2026-05-23] detrecord/irectotal retired — see state%atmosphere%{detrecord,irectotal}
      ! [GR-CROP-DVS] nmetdetail retired — see state%cfg%meteo%nmetdetail
      ! [GR-CROP C12] retired 2026-05-14 — moved to state%atmosphere%nmrain (Arc 8)
      ! integer   nmrain             ! Number of rain event records (-)
      integer   nofd               ! number of days for running average Tmin (-)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%rainrec (ADR 0041)
      ! integer   rainrec            ! moved to state%timecontrol%rainrec
      ! [GR-CROP-DVS] swdivide retired — see state%cfg%meteo%swdivide
      ! [GR-CROP-DVS] swetr retired — see state%cfg%meteo%swetr
      integer   swetsine           ! Switch: 0 = Tp and Ep uniform during a day; 1 = Tp and Ep are distributed as sine waves during a day
      ! [GR-CROPWS] swinter retired — see state%crop%common%swinter
      ! [GR-CROP-DVS] swmetdetail retired — see state%cfg%meteo%swmetdetail
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%swmeteo (ADR 0041)
      ! integer   swmeteo            ! moved to state%timecontrol%swmeteo
      ! [GR-CROP-DVS] swrain retired — see state%cfg%meteo%swrain
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
      ! [GR-CROP-DVS] alt retired — see state%cfg%meteo%alt
      ! [GR-CROP-DVS] altw retired — see state%cfg%meteo%altw
      ! [GR-ATM 2026-05-23] angstroma/b retired — see state%atmosphere%angstrom{a,b}
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
      ! [GR-ATM 2026-05-23] daylp retired — see state%atmosphere%daylp
      ! [GR-ATM 2026-05-23] dethum/detrad/detrain/dettav/dettime/detwind retired — see state%atmosphere%det*   
      ! [GR-ATM 2026-05-23] dtEventRain retired — see state%atmosphere%dtEventRain
      ! real(8) :: empreva   ! Reduced soil evaporation flux according to empirical functions (L/T)
      real(8)   epot(96)           ! In case of detailed weather input, calculated Epot of each weather record (L/T)
      ! [GR-ATM 2026-05-23] cfevappond retired — see state%atmosphere%cfevappond
      ! [GR-ATM 2026-05-23] finterception retired — see state%atmosphere%finterception
      ! real(8) :: fprecnosnow  ! Ratio rain (excl. snow and rain on snow) / gross rain flux (-)
      ! real(8) :: grai        ! Daily gross rain flux (L/T), without rain on snow
      ! real(8) :: graidt      ! Gross precipitation flux during iteration timesteps (L/T)
      real(8)   grain(96)          ! In case of detailed weather input, gross rain flux of each weather record (L/T)
      ! real(8) :: ievap       ! Intermediate amount of actual soil evaporation (L)
      ! real(8) :: inrai       ! Intermediate amount of net precipitation (L)
      ! real(8) :: ipeva       ! Intermediate amount of potential soil evaporation (L)
      ! real(8) :: iptra       ! Intermediate amount of potential transpiration (L)
      ! [GR-CROP-DVS] lat retired — see state%cfg%meteo%lat
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%metperiod (ADR 0041)
      ! real(8)   metperiod          ! moved to state%timecontrol%metperiod
      ! real(8) :: nraida      ! Daily average net precipitation flux (L/T)
      ! real(8) :: nraidt      ! Net precipitation flux during iteration timesteps (L/T)
      real(8)   nrain(96)          ! In case of detailed weather input, calculated netto rain of each weather record (L/T)
      ! real(8) :: peva        ! Potential soil evaporation flux (L/T)
      ! real(8) :: pevaday     ! Potential soil evaporation of one day (L)
      ! real(8) :: ptra        ! Potential transpiration flux (L/T)
      ! real(8) :: ptraday     ! Potential transpiration of one day (L)
      ! [GR-ATM 2026-05-23] rad retired — see state%atmosphere%rad
      ! [GR-CROP C12] retired 2026-05-14 — moved to state%atmosphere%rainamount (Arc 8)
      ! real(8)   rainamount(mrain)  ! Array with short duration rainfall amounts (L)
      ! [GR-CROP-DVS] raintab retired — see state%cfg%meteo%raintab
      ! [GR-CROP C12] retired 2026-05-14 — moved to state%atmosphere%rainfluxarray (Arc 8)
      ! real(8)   rainfluxarray(mrain) ! Array with short duration rainfall intensities (L/T)
      ! [GR-CROP C12] retired 2026-05-14 — moved to state%atmosphere%raintimearray (Arc 8)
      ! real(8)   raintimearray(mrain) ! Array with times (T) at which rainfall intensity changes
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%rh (GR-ATM arc)
      ! real(8)   rh                 ! Relative air humidity (-)
      ! [GR-ATM 2026-05-23] tav retired — see state%atmosphere%Tav
      ! [GR-ATM C8] retired 2026-05-14 — moved to state%atmosphere%tavd (GR-ATM arc)
      ! real(8)   tavd               ! Average air temperature during day time (oC)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%timjan1 (ADR 0041)
      ! real(8)   timjan1            ! moved to state%timecontrol%timjan1
      ! [GR-ATM 2026-05-23] tmn retired — see state%atmosphere%tmn
      real(8)   tmnr               ! Average of minimum air temperature during past 7 days (oC)
      ! [GR-ATM 2026-05-23] tmx retired — see state%atmosphere%tmx
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
      ! [GR-ATM 2026-05-23] tsunrise_atm/tsunset_atm retired — see state%atmosphere%tsun{rise,set}_atm  
      integer   nod10_cn           ! Node at -10cm for CN runoff method - from meteoday.f90 CNmethod
      integer   icn_atm            ! Current position in CN time table - from meteoday.f90 CNmethod
      ! [GR-CROP-DVS] z10_cn retired — see state%atmosphere%z10_cn
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
      ! [SS-GR-FINAL D1] swirg retired — IRG output deleted; 0 consumers
      ! [SS-GR-FINAL D1] irg retired — IRG file handle, never opened; 0 consumers
      integer   irrigevent         ! Switch: 0 = no irrigation; 1 = fixed irrigation event; 2 = scheduled irrigation event
      integer   irtype(mairg)      ! Type of fixed irrigation: 0 = sprinkling irrigation; 1 = surface irrigation
      integer   isua               ! Switch for type of irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer   isuas              ! Switch for type of scheduled irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer   nirri              ! Number of irrigation event
      ! [SS-GR-FINAL D1] phormc retired — 0 consumers
      ! [GR-CROPWS] schedule retired — see state%crop%common%schedule
      integer   swirfix            ! Switch for fixed irrigation: 0 = no applications prescribed; 1 = applications are prescribed
      integer   swcirrthres        ! Switch to allow over irrigation when a conc-threshold is exceeded: 0 = no; 1 = yes/allowed
      real(8)   cirrs              ! Solute concentration of irrigation water (M/L3)
      real(8)   cirrthres          ! Threshold value (M/L3) indicating the concentration that initiates over irrigation
      real(8)   dcrit              ! Depth (L) of sensor for soil water pressure head or water content
      real(8)   ditab(14)          ! Array with amount of under- or over-irrigation (L) as function of crop development stage
      real(8)   dwatab(14)         ! Array with maximum amounts of water depleted as function of crop development stage
      real(8)   fidtab(14)         ! Array with prescribed fixed irrigation depth (L) as function of crop development stage
      ! [GR-CROP-DVS] gird retired — see state%crop%gird
      real(8)   hcritab(14)        ! Array with minimum soil water pressure heads (L) as function of crop development stage
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater (ADR 0038)
      ! real(8)   igird              ! Intermediate depth of gross irrigation (L)
      ! real(8)   inird              ! Intermediate depth of net irrigation (L)
      real(8)   irconc(mairg)      ! Array with irrigation concentrations (M/L3) in case of fixed irrigation
      real(8)   irdate(mairg)      ! Array with fixed irrigation dates
      real(8)   irdepth(mairg)     ! Array with fixed irrigation depths (L)
      ! [GR-ATM 2026-05-23] nird retired — see state%atmosphere%nird
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
      ! [SS-GR-FINAL D1] FlIrrigationOutput retired — 0 consumers
      integer   dayfix             ! days since last irrigation event

! --- start of growth grass      
      ! [GR-CROPWS] swtsum retired — see state%crop%grass%swtsum
      integer   tsumtime           ! time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      real(8)   tsumtemp           ! temperature limit to initiate grass growth  [0.0..20.0 grC, R]
      real(8)   tsumdepth          ! depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
      
! --- crop variables
      ! [SS-GR-FINAL D1] swCrop retired — all callers use state%crop; 0 consumers
      logical   flCropCalendar     ! Flag indicating that crop season is active (but currently might be bare or cropped)
      logical   flCropEmergence    ! Flag indicating period from crop emergence until harvest
      logical   flCropHarvest      ! Flag indicating period from crop harvest until the end of crop season
      logical   flCropReadFile     ! Flag indicating reading of input.crp
      logical   flCropOpenFile     ! Flag indicating to create output.crp
      logical   flCropOutput       ! Flag indicating writing of output.crp

      ! [GR-ATM 2026-05-23] croptype retired — see state%crop%common%croptype
      integer   swcrp              ! Switch for output file *.CRP with daily crop output: 0 = no; 1 = yes
      integer   crp                ! Internal number of crop output file *.CRP
      integer   daycrop            ! Number of days that a crop exists
      integer   icrop              ! Current crop number
      ! [GR-CROPWS] daygrowth retired — see state%crop%grass%daygrowth
      ! [GR-CROPWS] daygrowthpot retired — see state%crop%grass%daygrowthpot
      ! [GR-CROPWS] idaysgraz retired — see state%crop%grass%idaysgraz
      ! [GR-CROPWS] idaysgrazpot retired — see state%crop%grass%idaysgrazpot
      ! [GR-CROPWS] idev retired — see state%crop%common%idev
      ! [GR-CROPWS] idregr retired — see state%crop%grass%idregr
      ! [GR-CROPWS] idregrpot retired — see state%crop%grass%idregrpot
      integer   idsl               ! Switch for crop development before anthesis: 0 = depends on temperature; 
                                   !   1 = depends on temperature and day length; 2 = depends on temperature, day length and vernalisation factor
      ! [GR-CROPWS] iharvest retired — see state%crop%grass%iharvest
      ! [GR-CROPWS] ilvold retired — see state%crop%common%ilvold
      ! [GR-CROPWS] ilvoldpot retired — see state%crop%common%ilvoldpot
      ! [GR-CROPWS] iseqgm retired — see state%crop%grass%iseqgm
      ! [GR-CROPWS] iseqgmpot retired — see state%crop%grass%iseqgmpot
      integer   noddrz             ! Compartment number at bottom root zone (-)
      ! [GR-CROPWS] seqgrazmow retired — see state%crop%grass%seqgrazmow
      ! [GR-CROPWS] seqgrazmowpot retired — see state%crop%grass%seqgrazmowpot
      ! [GR-CROPWS] swcf retired — see state%crop%swcf
      ! [GR-CROPWS] swdrought retired — see state%crop%common%swdrought
      ! [GR-CROPWS] swgc retired — see state%crop%common%swgc
      ! [SS-GR-FINAL D1] swjarvis retired — deprecated switch; 0 consumers
      ! [GR-CROPWS] swcompensate retired — see state%crop%common%swcompensate
      ! [GR-CROPWS] alphacrit retired — see state%crop%common%alphacrit
      ! [GR-CROPWS] dcritrtz retired — see state%crop%common%dcritrtz
      ! [GR-CROPWS] swstressor retired — see state%crop%common%swstressor
      ! [GR-CROPWS] swrootradius retired — see state%crop%common%swrootradius
      ! [GR-CROPWS] swsalinity retired — see state%crop%common%swsalinity
      ! [GR-ATM 2026-05-23] atmtr retired — see state%atmosphere%atmtr
      real(8)   agerm              ! Coefficient a  of germination
      real(8)   cgerm              ! Coefficient c  of germination
      real(8)   bgerm              ! Coefficient b  of germination
      ! [GR-CROPWS] adcrh retired — see state%crop%common%adcrh
      ! [GR-CROPWS] adcrl retired — see state%crop%common%adcrl
      ! [GR-CROPWS] air_filled_root_por retired — see state%crop%common%air_filled_root_por
      ! [GR-CROP-DVS] albedo retired — see state%crop%common%albedo
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%alpJvLier (ADR 0036)
      ! real(8)   alpJvLier          ! Uniform drought reduction factor based on concept Jong van Lier (-)
      ! [GR-CROPWS] amaxtb retired — see state%crop%common%amaxtb
      ! [GR-CROP-DVS] pgasspot retired — see state%crop%wofost%pgasspot
      ! [GR-CROP-DVS] pgass retired — see state%crop%wofost%pgass
      integer   swpotrelmf         ! Calculation of potential yield
      ! [GR-CROP-DVS] relmf retired — see state%crop%grass%relmf
      real(8)   avevaptb(2*magrs)  ! Gash interception model: average evaporation intensity during shower (-) as function of time (T)
      real(8)   avprectb(2*magrs)  ! Gash interception model: average rainfall intensity (-) as function of time (T)
      ! [GR-CROP-DVS] cf retired — see state%crop%common%cf
      ! [GR-CROPWS] cfeic retired — see state%crop%fixed%cfeic
      ! [GR-CROPWS] gctb retired — see state%crop%fixed%gctb
      ! [GR-CROPWS] cftb retired — see state%crop%fixed%cftb
      ! [GR-CROPWS] cfeictb retired — see state%crop%fixed%cfeictb
      ! [GR-CROP-DVS] ch retired — see state%crop%common%ch
      ! [GR-CROPWS] chtb retired — see state%crop%fixed%chtb
      real(8)   c_mroot            ! Maintenance coefficient of root [0.0..1.0 kg O2/kg/d, R]
      ! [GR-CROP-DVS] cofab retired — see state%crop%cofab
      real(8)   cropend(macrop)    ! Array with crop end dates
      real(8)   cropstart(macrop)  ! Array with crop start dates
      ! [GR-CROPWS] cumdens retired — see state%crop%common%cumdens
      ! [GR-CROPWS] cuptgraz retired — see state%crop%common%cuptgraz
      ! [GR-CROPWS] cuptgrazpot retired — see state%crop%common%cuptgrazpot
      ! [GR-CROPWS] cvl retired — see state%crop%common%cvl
      ! [GR-CROPWS] cvo retired — see state%crop%common%cvo
      ! [GR-CROPWS] cvr retired — see state%crop%common%cvr
      ! [GR-CROPWS] cvs retired — see state%crop%common%cvs
      ! [GR-CROP-DVS] cwdm retired — see state%crop%wofost%cwdm
      ! [GR-CROP-DVS] cwdmpot retired — see state%crop%wofost%cwdmpot
      ! [GR-ATM 2026-05-23] difpp retired — see state%atmosphere%difpp
      real(8)   dlc                ! Shortest day length (T) for any crop development
      real(8)   dlo                ! Minimum day length (T) for optimal crop development
      ! [GR-CROPWS] dry_mat_cont_roots retired — see state%crop%common%dry_mat_cont_roots
      ! [GR-ATM 2026-05-23] dsinbe retired — see state%atmosphere%dsinbe
      ! [GR-CROPWS] dtsmtb retired — see state%crop%common%dtsmtb
      ! [GR-CROP-DVS] dvs retired — see state%crop%common%dvs
      ! [GR-CROPWS] dvsend retired — see state%crop%common%dvsend
      ! [GR-CROP-DVS] dwlv retired — see state%crop%wofost%dwlv
      ! [GR-CROP-DVS] dwlvCrop retired — see state%crop%wofost%dwlvCrop
      ! [GR-CROP-DVS] dwlvSoil retired — see state%crop%wofost%dwlvSoil
      ! [GR-CROP-DVS] dwlvpot retired — see state%crop%wofost%dwlvpot
      ! [GR-CROP-DVS] dwrt retired — see state%crop%wofost%dwrt
      ! [GR-CROP-DVS] dwrtpot retired — see state%crop%wofost%dwrtpot
      ! [GR-CROP-DVS] dwso retired — see state%crop%wofost%dwso
      ! [GR-CROP-DVS] dwst retired — see state%crop%wofost%dwst
      ! [GR-CROP-DVS] dwstpot retired — see state%crop%wofost%dwstpot
      ! [GR-CROPWS] eff retired — see state%crop%common%eff
      real(8)   f_senes            ! Reduction factor for senescence, used for maintenance respiration [0..1.0 -, R]
      ! [GR-CROPWS] fltb retired — see state%crop%common%fltb
      ! [GR-CROPWS] fotb retired — see state%crop%common%fotb
      ! [GR-CROPWS] frtb retired — see state%crop%common%frtb
      ! [GR-CROPWS] fstb retired — see state%crop%common%fstb
      real(8)   gasst              ! Total gross assimilation for actual crop (kg/ha)
      real(8)   gasstpot           ! Total gross assimilation for potential crop (kg/ha)
      ! [GR-ATM 2026-05-23] gc retired — see state%crop%common%gc
      ! [GR-CROP-DVS] HarLosOrm_tot retired — see state%crop%common%HarLosOrm_tot
      real(8)   hdrygerm           ! Criterium Hdry of germination
      real(8)   hwetgerm           ! Criterium Hwet of germination
      ! [GR-CROPWS] hlim1 retired — see state%crop%common%hlim1
      ! [GR-CROPWS] hlim2l retired — see state%crop%common%hlim2l
      ! [GR-CROPWS] hlim2u retired — see state%crop%common%hlim2u
      ! [GR-CROPWS] hlim3h retired — see state%crop%common%hlim3h
      ! [GR-CROPWS] hlim3l retired — see state%crop%common%hlim3l
      ! [GR-CROPWS] hlim4 retired — see state%crop%common%hlim4
      ! [GR-CROP-DVS] kdif retired — see state%crop%kdif
      ! [GR-CROP-DVS] kdir retired — see state%crop%kdir
      ! [GR-CROP-DVS] lai retired — see state%crop%lai
      ! [GR-CROPWS] laiem retired — see state%crop%common%laiem
      ! [GR-CROPWS] laiexp retired — see state%crop%common%laiexp
      ! [GR-CROPWS] laiexppot retired — see state%crop%common%laiexppot
      ! [GR-CROPWS] glaiex retired — see state%crop%common%glaiex
      ! [GR-CROPWS] glaiexpot retired — see state%crop%common%glaiexpot
      ! [GR-CROPWS] laimax retired — see state%crop%common%laimax
      ! [GR-CROP-DVS] laipot retired — see state%crop%common%laipot
      ! [GR-CROPWS] lv retired — see state%crop%common%lv
      ! [GR-CROPWS] lvpot retired — see state%crop%common%lvpot
      ! [GR-CROPWS] lvage retired — see state%crop%common%lvage
      ! [GR-CROPWS] lvagepot retired — see state%crop%common%lvagepot
      real(8)   max_resp_factor    ! Ratio root total respiration / maintenance respiration [1..5.0 -, R]
      ! [GR-CROP-DVS] mowrest retired — see state%crop%grass%mowrest
      ! [GR-CROPWS] dewrest retired — see state%crop%grass%dewrest
      real(8)   mrest              ! Total maintenance respiration for actual crop (kg/ha)
      real(8)   mrestpot           ! Total maintenance respiration for potential crop (kg/ha)
      real(8)   mrftb(2*magrs)     ! Array with ratio root total respiration / maintenance respiration as function of DVS (kg/m3)
      ! [GR-CROPWS] perdl retired — see state%crop%common%perdl
      real(8)   pfreetb(2*magrs)   ! Gash interception model: free throughfall coefficient (-) as function of time (T)
      real(8)   pstemtb(2*magrs)   ! Gash interception model: stem flow coefficient (-) as function of time (T)
      ! [GR-ATM 2026-05-23] siccaptb retired — see state%crop%common%siccaptb
      ! [GR-CROP-DVS] fimin retired — see state%atmosphere%fimin
      ! [GR-CROP-DVS] siccapact retired — see state%atmosphere%siccapact
      ! real(8) :: sicact        ! amount of water stored on canopy (cm) — [SS-ATM] retired 2026-05-11
      real(8)   siccaplai          ! interception storage per unit of LAI (cm/LAI)
      ! [GR-CROPWS] q10 retired — see state%crop%common%q10
      real(8)   q10_microbial      ! Relative increase in microbial respiration at temperature increase of 10 �C [1.0..4.0 -, R]
      real(8)   q10_root           ! Relative increase in root respiration at temperature increase of 10 �C [1.0..4.0 -, R]
      ! [GR-CROPWS] rdrrtb retired — see state%crop%common%rdrrtb
      ! [GR-CROPWS] rdrstb retired — see state%crop%common%rdrstb
      real(8)   reltr              ! relative transpiration factor that reduces crop growth (-)
      ! [GR-CROPWS] rfsetb retired — see state%crop%common%rfsetb
      ! [GR-CROPWS] rgrlai retired — see state%crop%common%rgrlai
      real(8)   rid                ! Real day number of detailed grass crop (d)
      ! [GR-CROPWS] rml retired — see state%crop%common%rml
      ! [GR-CROPWS] rmo retired — see state%crop%common%rmo
      ! [GR-CROPWS] rmr retired — see state%crop%common%rmr
      ! [GR-CROPWS] rms retired — see state%crop%common%rms
      real(8)   rootcoefa          ! Defines relative distance at which mean soil water content occurs between roots
      real(8)   rooteff            ! Root system efficiency factor [0..1.0 -, R]
      real(8)   rootradius         ! Root radius drought stress (cm)
      ! [GR-CROPWS] root_radiusO2 retired — see state%crop%common%root_radiusO2
      ! [GR-CROP-DVS] rsc retired — see state%crop%common%rsc
      ! [GR-CROPWS] rsw retired — see state%crop%common%rsw
      real(8)   scanopytb(2*magrs) ! Gash interception model: storage capacity of canopy (-) as function of time (T)
      real(8)   shape_factor_rootr ! Shape factor for exponential decrease of root respiration rate with depth [0..1.0 -, R]
      ! [GR-CROPWS] sla retired — see state%crop%common%sla
      ! [GR-CROPWS] slapot retired — see state%crop%common%slapot
      ! [GR-CROPWS] slatb retired — see state%crop%common%slatb
      ! [GR-CROPWS] spa retired — see state%crop%common%spa
      ! [GR-CROPWS] span retired — see state%crop%common%span
      ! [GR-CROPWS] spec_weight_root_tissue retired — see state%crop%common%spec_weight_root_tissue
      real(8)   specific_resp_humus ! Respiration rate of humus at 25 �C [0.0..1.0 kg O2/kg C/d, R] 
      ! [GR-CROPWS] srl retired — see state%crop%common%srl
      ! [GR-CROPWS] ssa retired — see state%crop%common%ssa
      real(8)   tadw               ! Dry weight of plant minus roots of actual growth (kg/ha)
      real(8)   tadwpot            ! Dry weight of plant minus roots of potential growth (kg/ha)
      ! [GR-CROP-DVS] tagp retired — see state%crop%wofost%tagp
      ! [GR-CROP-DVS] tagppot retired — see state%crop%wofost%tagppot
      ! [GR-CROP-DVS] tagpt retired — see state%crop%wofost%tagpt
      ! [GR-CROP-DVS] tagptpot retired — see state%crop%wofost%tagptpot
      ! [GR-CROPWS] tbase retired — see state%crop%common%tbase
      ! [GR-CROPWS] tdwi retired — see state%crop%common%tdwi
      ! [GR-CROPWS] tmnftb retired — see state%crop%common%tmnftb
      ! [GR-CROPWS] tmpftb retired — see state%crop%common%tmpftb
      ! [GR-CROP-DVS] tsum retired — see state%crop%common%tsum
      ! [GR-CROPWS] tsumam retired — see state%crop%common%tsumam
      ! [GR-CROPWS] tsumea retired — see state%crop%common%tsumea
      real(8)   tsumemeopt         ! Temperature sum for crop emergence under optimal conditions
      real(8)   tsumgerm           ! Temperature sum during germination
      real(8)   TBASEM             ! Lower threshold temp. for emergence (C)
      real(8)   TEFFMX             ! max. eff. temp. for emergence (C)

      ! [GR-CROPWS] var_a retired — see state%crop%common%var_a
      real(8)   w_root_ss          ! Dry weight of roots at soil surface [0.0..10.0 kg/m3, R]
      real(8)   wiltpoint          ! Minimum pressure head at interface soil-root (cm)

      ! [GR-CROP-DVS] wlv retired — see state%crop%wofost%wlv
      ! [GR-CROP-DVS] wlvpot retired — see state%crop%wofost%wlvpot
      real(8)   wrtb(2*magrs)      ! Array with dry weight of root at soil surface as function of DVS (kg/m3)
      ! [GR-CROP-DVS] wrt retired — see state%crop%wofost%wrt
      ! [GR-CROP-DVS] wrtpot retired — see state%crop%wofost%wrtpot
      real(8)   wrtmin             ! Minimum dry weight of plant root at relative depth (1% of the initial value)
      real(8)   gwrt               ! Growth of dry weight of plant root (kg/ha)
      ! [GR-CROP-DVS] wso retired — see state%crop%wofost%wso
      ! [GR-CROP-DVS] wsopot retired — see state%crop%wofost%wsopot
      ! [GR-CROP-DVS] wst retired — see state%crop%wofost%wst
      ! [GR-CROP-DVS] wstpot retired — see state%crop%wofost%wstpot
      logical   flanthesis         ! Flag indicating anthesis stage of a crop
      ! [GR-CROPWS] flHarvest retired — see state%crop%grass%flHarvest
      logical   flHarvestDay        ! Flag indicating that current day is harvest day
      ! [GR-CROPWS] flHarvestpot retired — see state%crop%grass%flHarvestpot
      ! [GR-CROPWS] flGrazing retired — see state%crop%grass%flGrazing
      ! [GR-CROPWS] flGrazingpot retired — see state%crop%grass%flGrazingpot
      character(len=40) cropfil(macrop)   ! Array with names of crop files
      character(len=80) pathcrop          ! Path to folder with crop input files
      ! [GR-CROP-DVS] inifil retired — dead (no readers, no writers)

!     Harvest Grass
      ! [GR-CROPWS] dmmowtb retired — see state%crop%grass%dmmowtb
      ! [GR-CROPWS] dmgrztb retired — see state%crop%grass%dmgrztb
      ! [GR-CROPWS] dateharvest retired — see state%crop%grass%dateharvest
      ! [GR-CROPWS] DelayRegrowthTab retired — see state%crop%grass%DelayRegrowthTab
      ! [GR-CROPWS] lsda retired — see state%crop%grass%lsda
      ! [GR-CROPWS] DaysGrazingtab retired — see state%crop%grass%daysgrazingtab
      ! [GR-CROPWS] UptGrazingtab retired — see state%crop%grass%uptgrazingtab
      ! [GR-CROPWS] LossGrazingtab retired — see state%crop%grass%lossgrazingtab
      ! [GR-CROPWS] lossmowtab retired — unused
      ! [GR-CROPWS] lossgrztab retired — unused
      
      ! [GR-CROP-DVS] cropstartpot retired — see state%crop%grass%cropstartpot
      ! [GR-CROP-DVS] cropendpot retired — see state%crop%grass%cropendpot
      ! [GR-CROP-DVS] cropstartact retired — see state%crop%grass%cropstartact
      ! [GR-CROP-DVS] cropendact retired — see state%crop%grass%cropendact
      ! [GR-CROPWS] flhrvendpot retired — see state%crop%grass%flhrvendpot
      ! [GR-CROPWS] flhrvendact retired — see state%crop%grass%flhrvendact
      ! [GR-CROPWS] pmowdm retired — see state%crop%grass%pmowdm
      ! [GR-CROPWS] mowdm retired — see state%crop%grass%mowdm
      ! [GR-CROPWS] pgrzdm retired — see state%crop%grass%pgrzdm
      ! [GR-CROPWS] grzdm retired — see state%crop%grass%grzdm
      ! [GR-CROP-DVS] plossdm retired — see state%crop%wofost%plossdm
      ! [GR-CROP-DVS] lossdm retired — see state%crop%wofost%lossdm
      
! --- rooting      
      ! [GR-CROP-DVS] swrdc retired — see state%crop%common%swrdc
      ! [GR-CROPWS] rdctb retired — see state%crop%common%rdctb
      ! [GR-CROP-DVS] swrd retired — see state%crop%common%swrd
      ! [GR-CROPWS] rdtb retired — see state%crop%common%rdtb
      ! [GR-CROPWS] rlwtb retired — see state%crop%common%rlwtb
      ! [GR-CROPWS] wrtmax retired — see state%crop%common%wrtmax
      ! [GR-CROP-DVS] swdmi2rd retired — see state%crop%common%swdmi2rd
      ! [GR-CROP-DVS] rdi retired — see state%crop%common%rdi
      ! [GR-CROP-DVS] rri retired — see state%crop%common%rri
      ! [GR-CROP-DVS] rdc retired — see state%crop%common%rdc
      real(8)   rdmax              ! Maximum rooting depth in soil profile (L)
      ! [GR-CROP-DVS] rd retired — see state%crop%common%rd
      ! [GR-CROP-DVS] rdpot retired — see state%crop%common%rdpot
      ! [GR-CROP-DVS] rdm retired — see state%crop%common%rdm

! --- CO2
      ! [GR-ATM 2026-05-23] flCO2 retired — see state%atmosphere%flco2 (dormant, no TOML wiring yet)
      ! [GR-CROP-DVS] fco2amax retired — see state%crop%wofost%fco2amax
      ! [GR-CROP-DVS] fco2eff retired — see state%crop%wofost%fco2eff
      ! [GR-CROP-DVS] fco2tra retired — see state%crop%wofost%fco2tra
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
      ! [GR-CROP-DVS] dwbl retired — see state%crop%wofost%dwbl
      ! [GR-CROP-DVS] dwblpot retired — see state%crop%wofost%dwblpot
      real(8)   fbl                ! Dry weight fraction partitioned to flowers (-)
      ! [GR-CROPWS] fbltb retired — see state%crop%common%fbltb
      real(8)   pld
      real(8)   remoc
      ! [GR-CROP-DVS] plwt retired — see state%crop%wofost%plwt
      ! [GR-CROP-DVS] plwti retired — see state%crop%wofost%plwti
      ! [GR-CROP-DVS] wbl retired — see state%crop%wofost%wbl
      ! [GR-CROP-DVS] wblpot retired — see state%crop%wofost%wblpot

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
      ! [SS-GR-FINAL D1] amFERT retired — 0 consumers

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
      ! [GR-CROP-DVS] MaxBackTr retired — see state%cfg%simulation%numerical%MaxBackTr
      integer   Itnumb(100,2)      ! Iteration number statistics [soilhydraulics.f90, timecontrol_mod.f90]
      ! [GR-CROP-DVS] critdevh1cp retired — see state%cfg%simulation%numerical%critdevh1cp
      ! [GR-CROP-DVS] critdevh2cp retired — see state%cfg%simulation%numerical%critdevh2cp
      ! [GR-CROP-DVS] critdevponddt retired — see state%cfg%simulation%numerical%critdevponddt
      logical   fldumpconvcrit     ! flag to generate additional output about convergence-warnings from subr Headcalc
      logical   flwarn_hc          ! Headcalc warning flag (previously SAVE variable)
      integer   iwarn_hc           ! Headcalc warning counter (previously SAVE variable)
      ! [SS-GR-FINAL D1] nstep_hc retired — 0 consumers
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
      ! [SS-GR-FINAL D1] swcapriseoutput retired — 0 consumers

      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%H0max/k1max/q0 (GR-BH arc)
      ! real(8)   h0max
      ! real(8)   k1max
      ! real(8)   q0

      ! [SS-GR-FINAL D1] afo retired — AFO output deleted; 0 consumers
      ! [SS-GR-FINAL D1] aun retired — AUN output deleted; 0 consumers
      ! [SS-GR-FINAL D1] bal retired — BAL output deleted; 0 consumers
      ! [SS-GR-FINAL D1] blc retired — BLC output deleted; 0 consumers
      ! [SS-GR-FINAL D1] bma retired — BMA macropore balance never opened; 0 consumers
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%bpegwl (ADR 0038)
      ! integer   bpegwl             ! Node at bottom of perched groundwater
      integer   botcom(maho)       ! Array with number of bottom compartments in each soil layer
      integer   dra                ! Internal number of drainage input file *.DRA
      ! [GR-DRA 2026-05-23] dramet retired — see state%drainage%dramet
                                   !                              3 = drainage/infiltration resistance
      ! [SS-GR-CROPRT A3] swinc retired — always 0; outinc calls dropped (no config field)
      integer   inc                ! Internal number of output file *.INC with incremental water balance data
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%indeks (ADR 0038)
      ! integer   indeks(macp)       ! Index denoting wetting or drying curve in case of hysteresis: 1 = wetting; -1 = drying
      ! [GR-DRA 2026-05-23] ipos retired — see state%drainage%ipos
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
      ! [SS-GR-FINAL D1] swstr retired — STR output deleted; 0 consumers
      ! [SS-GR-FINAL D1] str retired — STR file handle; 0 consumers (swap_csv_output str is a local char var)
      integer   sw2                ! Switch for prescribed bottom flux: 1 = sine function; 2 = table
      integer   sw3                ! Switch for prescribed hydraulic head of deep aquifer: 1 = sine function; 2 = table
      integer   sw4                ! Switch for extra groundwater flux as function of time: 0 = no extra flux; 1 = include extra flux
      ! [SS-GR-CROPRT A3] swcsv retired — migrated to config%output_csv%enabled (readers use config now)
      character(len=1024) InList_csv   ! character string with comma-separated list of variables for CSV output
      ! [SS-GR-CROPRT A3] swcsv_tz retired — migrated to config%output_csv%enabled_tz
      character(len=1024) InList_csv_tz   ! character string with comma-separated list of variables for CSV output
      real(8)   tz_z1_z2(2)        ! Depth range for time-depth CSV output (default: top soil profile, bottom soil profile)
      ! [SS-GR-FINAL D1] swafo retired — AFO output deleted; 0 consumers
      ! [SS-GR-FINAL D1] swaun retired — AUN output deleted; 0 consumers
      ! [SS-GR-FINAL D1] swbal retired — BAL output deleted; 0 consumers
      ! [SS-GR-FINAL D1] swblc retired — BLC output deleted; 0 consumers
      ! [SS-GR-FINAL D1] swsba retired — SBA output deleted; 0 consumers
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%swbotb_runtime (GR-BH arc)
      ! integer   swbotb             ! Switch for bottom boundary condition (see *.SWP input file for overview)
      integer   swbotb3Impl        ! Switch for implicit solution with lower boundary option 3 (Cauchy): 0 = explicit, 1 = implicit
      ! [GR-BND 2026-05-23] SwBotb3ResVert retired — see state%soilwater%swbotb3resvert
      integer   swcfbs             ! Switch for use of coefficient CFBS to convert potential ET into potential E: 0 = no; 1 = yes
      integer   swdiscrvert        ! Switch to convert vertical discretization for soil water quality models: 0 = no; 1 = yes
      ! [GR-BH Task 37] swdivd retired — moved to state%drainage%swdivd
      ! [GR-BH Task 37] swdivdinf retired — moved to state%drainage%swdivdinf
      ! [GR-DRA 2026-05-23] swdislay retired — see state%drainage%swdislay
      ! [GR-DRA 2026-05-23] swtopdislay retired — see state%drainage%swtopdislay
      ! [GR-CROP-DVS] swdra retired — see state%surfacewater%swdra
                                   !                                            2 = simulate drainage and surface water
      integer   swfrost            ! Switch for reduction of hydraulic conductivity in case of frost: 0 = no; 1 = yes
      integer   swhyst             ! Switch for hysteresis of soil moisture retention function: 0 = no; 1 = yes
      ! [GR-SOL 2026-05-24] swinco retired — see state%soilwater%swinco
                                   !                                             3 = final pressure heads from previous simulation
      ! [GR-CROP-DVS] SWkmean retired — see state%cfg%simulation%numerical%swkmean
                                   !                                            3 = unweighted geometric mean,4 = weighted geometric mean
                                   !                                            5 = unweighted harmonic mean, 6 = weighted harmonic mean
      ! [GR-CROP-DVS] SwkImpl retired — see state%cfg%simulation%numerical%swkimpl
      ! [GR-DRA 2026-05-23] swliminf retired — see state%drainage%swliminf
      ! [SS-GR-FINAL D5] swoutputmodflow retired — 0 consumers; MODFLOW output deleted
      ! [GR-CROPWS] swoxygen retired — see state%crop%common%swoxygen
      ! [GR-CROPWS] swoxygentype retired — see state%crop%common%swoxygentype
      ! [GR-BND 2026-05-23] swpondmx retired — see state%surfacewater%swpondmx
      integer   swqhbot            ! Switch for flux-groundwater level relationship: 1 = exponential function; 2 = tabular function
      integer   swcofqhc           ! Switch for additional flux added to exponential flux-groundwater level relationship: 0 = no, 1 = yes
      ! [GR-ATM 2026-05-23] swredu retired — see state%atmosphere%swredu
                                   !                                           2 = use function of Boesten/Stroosnijder
      ! [GR-CROP-DVS] swsophy retired — see state%soilwater%swsophy
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      integer   ientrytablay(maho,0:matabentries) ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer

      integer   swtopsub           ! Switch for topsoil or subsoil: 1 = topsoil, 2 = subsoil
      
      ! [SS-GR-CROPRT A3] swrum retired — always 0 (no config field); outrume calls dropped
      ! [SS-GR-FINAL D1] swini retired — 0 consumers
      ! [SS-GR-CROPRT D1] swend retired — end-state dump switch; ADR 0009 zeroed permanently; state field deleted
      ! [SS-GR-FINAL D1] swwba retired — WBA output deleted; 0 consumers
      ! [SS-GR-FINAL D1] swvap retired — VAP output deleted; 0 consumers
      ! [SS-GR-FINAL D1] vap retired — VAP file handle; 0 consumers
      ! [SS-GR-FINAL D1] wba retired — WBA file handle; 0 consumers
      real(8)   aqamp              ! Amplitude of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(8)   aqave              ! Average hydraulic head in deep aquifer (L)
      real(8)   aqper              ! Period of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(8)   aqtmax             ! Time with maximum hydraulic head in deep aquifer (T)
      ! [SS-GR-FINAL D1] atop retired — reprofunctions array never read; 0 consumers
      ! [GR-DRA 2026-05-23] basegw retired — see state%drainage%basegw
      ! [GR-SOL 2026-05-24] bdens retired — see state%soilwater%bdens
      real(8)   c_top(macp)        ! Oxygen concentration at top of compartment(kg/m3)
      
      ! Oxygen stress per-node arrays (scaffolding from OxygenStress subroutine)
      real(8)   o2_d_soil_term1(macp)     ! Pre-calculated soil diffusion term1 per node
      real(8)   o2_d_soil_term2(macp)     ! Pre-calculated soil diffusion term2 per node
      real(8)   o2_gfp100(macp)           ! Gas-filled porosity * 100 per node
      real(8)   o2_capac_term(macp)       ! Water capacity term per node
      real(8)   o2_nmin1(macp)            ! N-1 per node for VG equation
      real(8)   o2_mplus1(macp)           ! M+1 per node for VG equation
      logical   o2_ini_stress             ! O2 stress initialization flag (initialized to .true. via data statement)
      
      ! [GR-CROP-DVS] cfbs retired — see state%crop%cfbs
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
      ! [GR-ATM 2026-05-23] cofred retired — see state%atmosphere%cofred
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
      ! [GR-DRA 2026-05-23] entres retired — see state%drainage%entres
      ! [GR-CROP-DVS] es0 retired — see state%crop%es0
      ! [GR-CROP-DVS] et0 retired — see state%crop%et0
      ! [GR-CROP-DVS] ew0 retired — see state%crop%ew0
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%evp (ADR 0038)
      ! real(8)   evp(macp)          ! Internal evaporation flux of top soil compartments (L/T)
      ! [GR-BH Task 37] FacDpthInf retired — moved to state%drainage%FacDpthInf
      ! [GR-DRA 2026-05-23] ftopdislay retired — see state%drainage%ftopdislay
      ! [GR-DRA 2026-05-23] geofac retired — see state%drainage%geofac
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%gwl (ADR 0038)
      ! real(8)   gwl                ! Groundwater level (L)
      ! [GR-CROP-DVS] gwlconv retired — see state%cfg%simulation%numerical%gwlconv
      real(8)   gwli               ! Groundwater level (L) at start of simulation
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   gwlinp             ! Prescribed groundwater level (L) for current time
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%gwlm1 (ADR 0038)
      ! real(8)   gwlm1              ! Groundwater level (L) at former time level
      ! [GR-BND 2026-05-23] gwltab retired — see state%soilwater%gwltab
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%h (ADR 0038)
      ! real(8)   h(macp)            ! Soil water pressure head (L)
      real(8)   h_enpr(macp)       ! Soil water Entry Pressure head for Modified MualemVanGenuchten curve (L)
      ! [GR-BND 2026-05-23] haqtab retired — see state%soilwater%haqtab
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%hatm (ADR 0038)
      ! real(8)   hatm               ! Pressure head of air (L) near the soil surface
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   hbot               ! Soil water pressure head (L) at bottom of soil column
      ! [GR-BND 2026-05-23] hbotab retired — see state%soilwater%hbotab
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
      ! [GR-SOL 2026-05-24] inpola/inpolb retired — see state%mesh%{inpola,inpolb}
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
      ! [SS-SWC] iqtdo/iqtup retired — iqinfmax stays (used in swap_csv_output); qinfmax retired
      real(8)   iqinfmax
      ! [SS-GR-FINAL D1] qinfmax retired — 0 consumers
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
      ! [GR-ATM 2026-05-23] ISsnowBeg retired — see state%atmosphere%ISsnowBeg
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%IThetaBeg (ADR 0038)
      ! real(8)   ithetabeg(macp)    ! Array with volumetric soil water contents (-) for each compartment at start of intermediate period
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%k (ADR 0038)
      ! real(8)   k(macp+1)          ! Array with soil hydraulic conductivity (L/T) for each numerical compartment
      ! [GR-DRA 2026-05-23] khbot retired — see state%drainage%khbot
      ! [GR-DRA 2026-05-23] khtop retired — see state%drainage%khtop
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%kmean (ADR 0038)
      ! real(8)   kmean(macp+1)      ! Array with mean soil hydraulic conductivity (L/T) at the interface of current and upper compartment
      real(8)   Kroot              ! Hydraulic, radial conductivity of root tissue (L/T)
      ! [GR-BH Task 36] retired 2026-05-13 — moved to state%soilwater%ksatfit/ksatexm (GR-BH arc)
      ! real(8)   ksatfit(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: fitted on VG based on lab data
      ! real(8)   ksatexm(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: examined in lab or field
      real(8)   ksatthr(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: to interpolate VG and Ksatexm
      real(8)   kstem              ! Conductance in the path from leaf to root xylem (/d)
      ! [GR-DRA 2026-05-23] kvbot retired — see state%drainage%kvbot
      ! [GR-DRA 2026-05-23] kvtop retired — see state%drainage%kvtop
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
      ! [GR-CROP-DVS] pondmx retired — see state%surfacewater%pondmx
      ! [GR-BND 2026-05-23] pondmxtab retired — see state%surfacewater%pondmxtab
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%q (ADR 0038)
      ! real(8)   q(macp+1)          ! Soil water flux between current compartment and upper compartment (L/T)
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   qbot               ! Water flux through bottom of simulated soil column (L/T)
      ! [GR-BND 2026-05-23] qbotab retired — see state%soilwater%qbotab
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! real(8)   qbot_nonfrozen     ! Water flux through bottom of non-frozen soil column (L/T)
      ! real(8)   qdra(Madr,macp)    ! Moved to drainage_state_t%qdra (ADR 0031)
      ! real(8)   qdrain(Madr)       ! Moved to drainage_state_t%qdrain (ADR 0031)
      real(8)   qdraincomp(macp)   ! Total lateral drainage flux (L/T) for each compartment
      ! [GR-DRA 2026-05-23] qdrtab retired — see state%drainage%qdrtab
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
      ! [SS-GR-FINAL D1] ResultsOxygenStress retired — 0 consumers 
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
      ! [GR-ATM 2026-05-23] rsigni retired — see state%atmosphere%rsigni
      ! [GR-ATM 2026-05-23] rsoil retired — see state%atmosphere%rsoil
      ! [GR-CROP-DVS] rsro retired — see state%surfacewater%rsro
      ! [GR-CROP-DVS] rsroexp retired — see state%surfacewater%rsroexp
      ! [GR-ATM 2026-05-23] swuseCN retired — see state%atmosphere%swusecn (dormant, no TOML wiring yet)
      ! [GR-CROP-DVS] wc_cor/CNrefTAB/CNtimTAB/iCNtab retired — see state%atmosphere
      ! [SS-SWC] retired 2026-05-12 — crunoffCN/irunoCN moved to state%soilwater (ADR 0038)
      ! [GR-ATM] retired 2026-05-23 — CNdry/CNwet/ThetaRef/wc10/Runoff_CN moved to state%atmosphere
      real(8)   Rxylem             ! Mean radius of xylem tube inside roots (L)
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%runon (ADR 0038)
      ! real(8)   runon              ! Water runon flux (L/T)
      ! [GR-BND 2026-05-23] runonarr retired — see state%soilwater%runonarr
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
      ! [GR-CROP-DVS] taccur retired — see state%cfg%simulation%numerical%taccur
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
      ! [GR-DRA 2026-05-23] zintf retired — see state%drainage%zintf
      ! real(8)   ztopdislay(Madr)   ! Moved to drainage_state_t%ztopdislay (ADR 0031)
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flDrain (ADR 0041)
      ! logical   fldrain            ! moved to state%timecontrol%flDrain
      logical   FlHydrLift         ! Flag indicating release of water from root to soil is allowed
      ! [SS-SWC] retired 2026-05-12 — moved to state%soilwater%fllowgwl (ADR 0038)
      ! logical   fllowgwl           ! Flag indicating precribed groundwater level below bottom soil column
      ! [GR-BND 2026-05-23] flrunon retired — see state%soilwater%flrunon
      ! [SS-BND] retired 2026-05-11 — boundary subsystem migrated to state%soilwater (ADR 0035)
      ! logical   ftoph              ! Flag indicating that the pressure head is prescribed at the soil surface
      character(len=16) drfil      ! Name of drainage input file
      character(len=80) pathdrain  ! Path to folder with drainage input files

! --- heat variables
      ! [SS-GR-FINAL D1] nheat retired — heat_state.f90 declares its own nheat; 0 consumers here
      integer   swbotbhea          ! Switch for bottom boundary condition: 1 = heat flux is zero; 2 = prescribed temperature
      integer   swtopbhea          ! Switch for top boundary condition: 1 = use air temperatures; 2 = read measured surface temperatures
      integer   swcalt             ! Switch for method of soil water heat flow simulation: 1 = analytical method; 2 = numerical method
      integer   swhea              ! Switch for simulation of soil heat flow: 0 = no; 1 = yes
      ! [SS-GR-CROPRT A3] swtem retired — always 0 (no config field); outtem calls dropped
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
      ! [GR-CROP-DVS] TePrRain retired — see state%atmosphere%TePrRain
      ! [GR-CROP-DVS] TePrSnow retired — see state%atmosphere%TePrSnow
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flSnow (ADR 0041)
      ! logical   flsnow             ! moved to state%timecontrol%flSnow

! --- solute variables
      ! [GR-SOL 2026-05-24] nconc retired — see state%solute%nconc
      ! [SS-GR-FINAL D1] sba retired — SBA output deleted; 0 consumers
      ! [GR-SOL 2026-05-24] swbr retired — see state%solute%swbr
      ! [GR-SOL 2026-05-24] swbotbc retired — see state%solute%swbotbc
      integer   swsolu             ! Switch for simulation of solute transport: 0 = no; 1 = yes
      ! [SS-GR-FINAL D3] swsp retired — 0 consumers; sorption switch never read outside init/c2v
      ! [SS-GR-CROPRT A1] AgeGwl1m retired — ADR 0032 (AgeTracer dead-code)
      ! [GR-SOL 2026-05-24] bexp retired — see state%solute%bexp
      ! [GR-SOL 2026-05-24] cdrain retired — see state%solute%cdrain (AgeTracer dead-body deleted)
      ! [GR-SOL 2026-05-24] cirr retired — see state%solute%cirr
      ! [GR-SOL 2026-05-24] cml retired — see state%solute%cml_init / state%solute%cml
      real(8)   cmsy(macp)         ! Array with dissolved + adsorbed solute concentration (M/L3 soil volume) in mobile region
! real(8)   cpond              ! Moved to solute_state_t%cpond (ADR 0032)
      ! [GR-SOL 2026-05-24] cpre retired — see state%solute%cpre
      ! [GR-SOL 2026-05-24] cref retired — see state%solute%cref
! real(8)   cseep              ! Moved to solute_state_t%cseep (ADR 0032)
      ! [GR-SOL 2026-05-24] cseeptab retired — see state%solute%cseeptab
! real(8)   csurf              ! Moved to solute_state_t%csurf (ADR 0032)
      ! [GR-SOL 2026-05-24] daquif retired — see state%solute%daquif
      ! [GR-SOL 2026-05-24] ddif retired — see state%solute%ddif
      ! [GR-SOL 2026-05-24] decpot retired — see state%solute%decpot
      ! [GR-SOL 2026-05-24] decsat retired — see state%solute%decsat
! real(8)   dectot             ! Moved to solute_state_t%dectot (ADR 0032)
! real(8)   imdectot           ! Moved to solute_state_t%imdectot (ADR 0032)
      ! [GR-SOL 2026-05-24] dtsolu retired — see state%solute%dtsolu (AgeTracer dead-body deleted)
      ! [GR-SOL 2026-05-24] fdepth retired — see state%solute%fdepth
      ! [GR-SOL 2026-05-24] frexp retired — see state%solute%frexp
      ! [GR-SOL 2026-05-24] gampar retired — see state%solute%gampar
      ! [SS-GR-CROPRT A1] icAgeBot retired — ADR 0032 (AgeTracer dead-code)
      ! [SS-GR-CROPRT A1] icAgeDra retired — ADR 0032 (AgeTracer dead-code)
      ! [SS-GR-CROPRT A1] icAgeRot retired — ADR 0032 (AgeTracer dead-code)
      ! [SS-GR-CROPRT A1] icAgeSur retired — ADR 0032 (AgeTracer dead-code)
      ! [GR-SOL 2026-05-24] isqbot retired — see state%solute%isqbot (AgeTracer dead-body deleted)
      ! [GR-SOL 2026-05-24] isqtop retired — see state%solute%isqtop (AgeTracer dead-body deleted)
      ! [GR-SOL 2026-05-24] kf retired — see state%solute%kf
      ! [GR-SOL 2026-05-24] kfsat retired — see state%solute%kfsat
      ! [GR-SOL 2026-05-24] ldis retired — see state%solute%ldis
      ! [GR-SOL 2026-05-24] poros retired — see state%solute%poros
      ! [GR-SOL 2026-05-24] rottot retired — see state%solute%rottot (AgeTracer dead-body deleted)
! real(8)   imrottot           ! Moved to solute_state_t%imrottot (ADR 0032)
      ! [GR-SOL 2026-05-24] rtheta retired — see state%solute%rtheta
      ! [GR-CROPWS] salthead retired — see state%crop%common%salthead
      ! [GR-CROPWS] saltmax retired — see state%crop%common%saltmax
      ! [GR-CROPWS] saltslope retired — see state%crop%common%saltslope
! real(8)   samcra             ! Moved to solute_state_t%samcra (ADR 0032)
      ! [GR-SOL 2026-05-24] samini retired — see state%solute%samini (AgeTracer dead-body deleted)
! real(8)   sampro             ! Moved to solute_state_t%sampro (ADR 0032)
! real(8)   solbal             ! Moved to solute_state_t%solbal (ADR 0032)
! real(8)   sqbot              ! Moved to solute_state_t%sqbot (ADR 0032)
! real(8)   imsqbot            ! Moved to solute_state_t%imsqbot (ADR 0032)
      ! [GR-SOL 2026-05-24] sqdra retired — see state%solute%sqdra (AgeTracer dead-body deleted)
! real(8)   imsqdra            ! Moved to solute_state_t%imsqdra (ADR 0032)
! real(8)   sqirrig            ! Moved to solute_state_t%sqirrig (ADR 0032)
! real(8)   imsqirrig          ! Moved to solute_state_t%imsqirrig (ADR 0032)
! real(8)   sqprec             ! Moved to solute_state_t%sqprec (ADR 0032)
! real(8)   imsqprec           ! Moved to solute_state_t%imsqprec (ADR 0032)
! real(8)   sqrap              ! Moved to solute_state_t%sqrap (ADR 0032)
! real(8)   sqsur              ! Moved to solute_state_t%sqsur (ADR 0032)
      ! [GR-SOL 2026-05-24] tscf retired — see state%solute%tscf
      ! [GR-SOL 2026-05-24] zc retired — see state%solute%zc_init
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%flSolute (ADR 0041)
      ! logical   flsolute           ! moved to state%timecontrol%flSolute
      ! [SS-GR-CROPRT A1] flAgeTracer retired — ADR 0032 (AgeTracer dead-code; always .false.)

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
      ! [SS-GR-FINAL D1] SwSoilShr retired — macropore shrinkage; 0 consumers
      ! [SS-GR-FINAL D1] ThetCrMp retired — macropore shrinkage; 0 consumers
      real(8) Z_Tp                 ! [retired-zero] kept: ArMpTp/ArMpSs gating
      real(8) CritUndSatVol        ! [retired-zero] kept: waterbalance watertable() arg
! --- macropore variables (selected — most retired by deletion)
      real(8) ArMpTp               ! [retired-zero] kept: ArMpSs assignment
      real(8) cQMpLatSs            ! [retired-zero] kept: soilhydraulics zero-write
      ! [GR-CROP-DVS] cQMpOutDrRap retired — dead
      ! [SS-GR-FINAL D1] dFdhMp retired — always 0; 0 consumers
      ! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%dtold (ADR 0041)
      ! real(8) dtold                ! moved to state%timecontrol%dtold
      real(8) DiPoCp(MaCp)         ! [retired-zero] kept: soilgrid refinement
      real(8) iQMpOutDrRap         ! [retired-zero] kept: swap_csv/swapoutput DRAINAGE accumulator
      ! [SS-GR-CROPRT A2] iQInTopLatDm1 retired — ADR 0040 (0 consumers; macropore BMA dead)
      ! [SS-GR-CROPRT A2] iQInTopLatDm2 retired — ADR 0040 (0 consumers; macropore BMA dead)
      ! [SS-GR-CROPRT A2] iQInTopVrtDm1 retired — ADR 0040 (0 consumers; macropore BMA dead)
      ! [SS-GR-CROPRT A2] iQInTopVrtDm2 retired — ADR 0040 (0 consumers; macropore BMA dead)
      real(8) IAvFrMpWlWtDm1(MaCp) ! [retired-zero] kept: soilgrid refinement
      real(8) IAvFrMpWlWtDm2(MaCp) ! [retired-zero] kept: soilgrid refinement
      real(8) iQExcMtxDm1Cp(MaCp)  ! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(8) iQExcMtxDm2Cp(MaCp)  ! [retired-zero] kept: soilgrid.f90 macropore redistribution
      real(8) iQOutDrRapCp(MaCp)   ! [retired-zero] kept: soilgrid.f90 macropore redistribution
      ! [SS-GR-CROPRT A2] IWaSrDm1Beg retired — ADR 0040 (macropore BMA writer; checkmassbal branch dropped)
      ! [SS-GR-CROPRT A2] IWaSrDm2Beg retired — ADR 0040 (macropore BMA writer; checkmassbal branch dropped)
      ! [SS-GR-CROPRT A2] WaSrDm1 retired — ADR 0040 (macropore BMA writer; checkmassbal branch dropped)
      ! [SS-GR-CROPRT A2] WaSrDm2 retired — ADR 0040 (macropore BMA writer; checkmassbal branch dropped)
      ! [SS-GR-CROPRT A2] WaSrDm1Ini retired — ADR 0040 (macropore wbalance term; branch dropped)
      ! [SS-GR-CROPRT A2] WaSrDm2Ini retired — ADR 0040 (macropore wbalance term; branch dropped)
      real(8) VlMpStDm1(MaCp)      ! [retired-zero] kept: soilgrid refinement
      real(8) VlMpStDm2(MaCp)      ! [retired-zero] kept: soilgrid refinement
      ! [SS-GR-CROPRT A2] IcTopMP retired — ADR 0040 (waterbalance macropore branch dropped)
      ! [SS-GR-FINAL D1] IDecMpRat retired — 0 consumers
      real(8) QExcMpMtx(MaCp)      ! [retired-zero] kept: waterbalance use clause
      real(8) QMaPo                ! [retired-zero] kept: waterbalance qbot term
      ! [GR-DRA 2026-05-23] QRapDra retired — see state%drainage%QRapDra
      ! [MACRO-RETIRE 2026-05-12] note: NumLevRapDra/RapDraReaExp/RapDraResRef
      ! belong to the drainage subsystem (set by drainage_config), NOT
      ! macropore — they were grouped here legacy-style. Keep them.
      ! [GR-DRA 2026-05-23] NumLevRapDra retired — see state%drainage%NumLevRapDra
      ! [SS-GR-FINAL D1] RapDraReaExp retired — only c2v writes; 0 consumers
      ! [SS-GR-FINAL D1] RapDraResRef retired — only c2v writes; 0 consumers
      logical FlDecMpRat           ! [retired-zero] kept: soilhydraulics convergence sentinel
      ! [SS-GR-CROPRT A2] flmacropore retired — ADR 0040 (always .false.; all guarded branches dropped)
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
      ! [SS-GR-FINAL D1] swswb retired — 0 consumers
      ! [SS-GR-FINAL D1] swdrf retired — 0 consumers
      ! [GR-CROP-DVS] swsrf retired — see state%cfg%surface_water%swsrf
      ! [GR-DRA 2026-05-23] swallo/swdtyp retired — see state%drainage%{swallo,swdtyp}
      ! [GR-BH Task 37] swnrsrf retired — moved to state%drainage%swnrsrf
      ! [GR-CROP-DVS] swsec retired — see state%cfg%surface_water%swsec
      ! [GR-DRA 2026-05-23] swqhr retired — see state%surfacewater%swqhr
      ! [GR-DRA 2026-05-23] nrpri retired — see state%surfacewater%nrpri
      ! [SS-GR-FINAL D1] nrsec retired — 0 consumers
      ! [GR-DRA 2026-05-23] nmper retired — see state%surfacewater%nmper
      ! [GR-DRA 2026-05-23] swman retired — see state%surfacewater%swman
      ! [GR-BH Task 37] SwTopnrsrf retired — moved to state%drainage%swtopnrsrf
      ! SS-SWST Phase 2 Task 11 C2: numadj/imper removed — state%surfacewater owns them.
      ! [SS-GR-FINAL D1] nqh retired — 0 consumers
      ! [SS-GR-FINAL D1] drf retired — 0 consumers
      ! [SS-GR-FINAL D1] swb retired — 0 consumers
      ! [GR-DRA 2026-05-23] nphase/nodhd retired — see state%surfacewater%{nphase,nodhd}
      ! numadj removed (surfacewater_state_t%numadj)
      ! [GR-DRA 2026-05-23] intwl retired — see state%surfacewater%intwl
      ! imper removed (surfacewater_state_t%imper)
      ! [GR-DRA 2026-05-23] nowltab retired — see state%drainage%nowltab
      real(8) widthr(Madr),taludr(Madr),rdrain(Madr),rsurfdeep
      real(8) rsurfshallow,rinfi(Madr),rentry(Madr),rexit(Madr)
      real(8) gwlinf(Madr)
      ! [GR-DRA 2026-05-23] wlptab retired — see state%surfacewater%wlptab
      ! [GR-DRA 2026-05-23] impend retired — see state%surfacewater%impend
      ! [GR-DRA 2026-05-23] wldip/wscap retired — see state%surfacewater%{wldip,wscap}
      ! [GR-DRA 2026-05-23] hbweir retired — see state%surfacewater%hbweir
      ! SS-SWST Phase 2 Task 11 C2: wlstar removed — state%surfacewater%wlstar owns it.
      ! [GR-DRA 2026-05-23] osswlm retired — see state%surfacewater%osswlm
      ! [GR-DRA 2026-05-23] alphaw/betaw retired — see state%surfacewater%{alphaw,betaw}
      real(8) wlp
      ! wlstar removed (surfacewater_state_t%wlstar)
      ! [SS-GR-FINAL D6] wls1_init retired — c2v write tombstoned in D4/C4; state%surfacewater owns it; 0 consumers
      ! [GR-DRA 2026-05-23] dropr retired — see state%surfacewater%dropr
      real(8) hdepth(mamp*mamte)
      ! [GR-DRA 2026-05-23] gwlcrit/hcrit/vcrit retired — see state%surfacewater%X
      real(8) hqhtab(mamp,mamte)
      ! [GR-DRA 2026-05-23] qqhtab retired — see state%surfacewater%qqhtab
      ! [GR-DRA 2026-05-23] wlsman retired — see state%surfacewater%wlsman
      ! SS-SWST Phase 2 Task 11 C2: sttab removed — state%surfacewater%sttab owns it.
      ! [GR-DRA 2026-05-23] wlstab retired — see state%surfacewater%wlstab
      ! sttab(22,2) removed (surfacewater_state_t%sttab)
      ! SS-SWST Phase 2 Task 11 C2: swstini/swst/wlsbak removed — state%surfacewater owns them.
      ! real(8) swstini   ! Moved to surfacewater_state_t%swstini
      ! real(8) swst      ! Moved to surfacewater_state_t%swst
      ! real(8) wlsbak(4) ! Moved to surfacewater_state_t%wlsbak
      ! SS-SWST Phase 2 Task 11 C2: cqdrd/cwsupp/cwout/wls removed — state%surfacewater owns them.
      ! [GR-DRA 2026-05-23] cofintfl/expintfl retired — see state%drainage%{cofintfl,expintfl}
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
      ! [SS-GR-FINAL D1] DayGerm retired — 0 consumers

      ! Harvest of crop growth
      ! [GR-CROPWS] swharv retired — see state%crop%common%swharv
      ! [SS-GR-FINAL D1] hmow retired — 0 consumers
      ! [SS-GR-FINAL D1] hgrz retired — 0 consumers
      ! [GR-CROPWS] zmow retired — unused
      ! [GR-CROPWS] zgrz retired — unused

! --- root development of dynamic crop growth and oxygen stress
      ! [GR-CROPWS] swWrtNonox retired — see state%crop%common%swWrtNonox
      ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%flWrtNonox (ADR 0036)
      ! logical   flWrtNonox         ! Flag indicating whether root development is retatarded by oxygen stress 
      ! [GR-CROPWS] aeratecrit retired — see state%crop%common%aeratecrit

      ! Initialize o2_ini_stress to .true. (needed for first call to OxygenStress)
      data o2_ini_stress /.true./

      end module variables
