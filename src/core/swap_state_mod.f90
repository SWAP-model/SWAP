! ==============================================================================
! SWAP State Module
! ==============================================================================
! This module defines the explicit state types for the SWAP model to enable:
!   - Multi-instance execution (multiple SWAP models in one process)
!   - Thread-safe parallelization  
!   - Future GPU offloading capability
!   - BMI-compatible state management
!
! Architecture:
!   swap_state_t          - Top-level container for all model state
!   io_handles_t          - File I/O handles (separate from simulation state)
!   *_state_t             - Domain-specific state types
!
! Usage:
!   type(swap_state_t) :: state
!   type(io_handles_t) :: io
!   call swap_state_init(state)
!   call io_handles_init(io)
!
! Author: SWAP Development Team
! Date: 2026-01-31
! ==============================================================================

module swap_state_mod
    use swap_log, only: log_info, to_str
    use soil_state_mod, only: soil_state_t, soil_state_init, soil_state_finalize
    use atmosphere_state_mod, only: atmosphere_state_t, atmosphere_state_init, atmosphere_state_finalize
    use boundary_state_mod, only: boundary_state_t, boundary_state_init, boundary_state_finalize
    use drainage_state_mod, only: drainage_state_t, drainage_state_init, drainage_state_finalize
    use surfacewater_state_mod, only: surfacewater_state_t, surfacewater_state_init, surfacewater_state_finalize
    use solute_state_mod, only: solute_state_t, solute_state_init, solute_state_finalize
    use heat_state_mod, only: heat_state_t, heat_state_init, heat_state_finalize
    use macropore_state_mod, only: macropore_state_t, macropore_state_init, macropore_state_finalize
    implicit none
    private

    ! Include array dimension parameters
    include 'arrays.fi'

    ! ===========================================================================
    ! Public types and procedures
    ! ===========================================================================
    public :: swap_state_t
    public :: io_handles_t
    public :: time_state_t
    public :: soil_state_t
    public :: atmosphere_state_t
    public :: crop_state_t
    public :: drainage_state_t
    public :: boundary_state_t
    public :: macropore_state_t
    public :: solute_state_t
    public :: heat_state_t
    public :: snow_state_t
    public :: irrigation_state_t
    public :: tillage_state_t
    public :: surfacewater_state_t
    public :: wofost_soil_state_t
    public :: oxygenstress_state_t
    
    ! Initialization procedures
    public :: swap_state_init
    public :: io_handles_init
    public :: drainage_state_init
    public :: surfacewater_state_init
    public :: macropore_state_init
    public :: irrigation_state_init
    public :: tillage_state_init
    public :: wofost_soil_state_init
    public :: oxygenstress_state_init
    
    ! Finalization procedures
    public :: swap_state_finalize
    public :: drainage_state_finalize
    public :: drain_state_finalize
    public :: surfacewater_state_finalize
    public :: macropore_state_finalize
    public :: irrigation_state_finalize
    public :: tillage_state_finalize
    public :: wofost_soil_state_finalize
    public :: oxygenstress_state_finalize

    ! ===========================================================================
    ! Time and Control State
    ! ===========================================================================
    type :: time_state_t
        ! Time stepping
        real(8) :: dt = 0.0d0              ! Current time step (d)
        real(8) :: dtmax = 0.0d0           ! Maximum time step (d)
        real(8) :: dtmin = 0.0d0           ! Minimum time step (d)
        real(8) :: t = 0.0d0               ! Time since start of calendar year (d)
        real(8) :: t1900 = 0.0d0           ! Time since 1900 (d)
        real(8) :: tcum = 0.0d0            ! Time since start of simulation (d)
        real(8) :: tend = 0.0d0            ! End date of simulation
        real(8) :: tstart = 0.0d0          ! Start date of simulation
        
        ! Counters
        integer :: daycum = 0              ! Day number from start of simulation
        integer :: daynr = 0               ! Day number of calendar year
        integer :: imonth = 0              ! Month number
        integer :: iyear = 0               ! Year number
        integer :: iyearm1 = 0             ! Previous year number
        integer :: isteps = 0              ! Number of time steps from start of day
        
        ! Output control
        integer :: ioutdat = 0             ! Counter of output date
        integer :: ioutdatint = 0          ! Counter of intermediate output date
        integer :: cntper = 0              ! Day number of intermediate period
        integer :: nprintday = 0           ! Number of output times during one day
        integer :: nprintcount = 0         ! Counter for output during a day
        integer :: period = 0              ! Length of prescribed output interval

        ! Output schedule configuration (loaded from .swp)
        integer :: swheader = 0            ! Header printing: 0=no, 1=yes
        integer :: swodat = 0              ! Extra output dates: 0=no, 1=yes
        integer :: swres = 0               ! Reset interval counter at year start
        integer :: swmonth = 0             ! Monthly intermediate output mode
        integer :: swscre = 0              ! Screen output mode
        real(8) :: outper = 0.0d0          ! Length of actual output interval
        real(8) :: outdat(maout) = 0.0d0   ! Output dates for balances
        real(8) :: outdatint(maout) = 0.0d0 ! Intermediate output dates
        
        ! Flags
        logical :: fldayend = .false.      ! End of day
        logical :: fldaystart = .false.    ! First time step of day
        logical :: flrunend = .false.      ! End of run
        logical :: flyearstart = .false.   ! Beginning of new year
        logical :: floutput = .false.      ! Time for output
        logical :: flbaloutput = .false.   ! Time for balance output
        logical :: flheader = .false.      ! Print header in output
        logical :: flzerocumu = .false.    ! Reset cumulative fluxes
        logical :: flzerointr = .false.    ! Reset intermediate fluxes
        logical :: fldecdt = .false.       ! Decrease time step
        logical :: fldecdtmin = .false.    ! Reset to minimum time step
        logical :: fldtmin = .false.       ! Time step equals minimum
        logical :: fldtreduce = .false.    ! Time step reduction flag
        logical :: flprintdt = .false.     ! Output every dt
        logical :: flprintshort = .false.  ! Several output times during day
        logical :: floutputshort = .false. ! Time for output during day
        
        ! Shared simulation flag
        logical :: flSwapShared = .false.  ! Simultaneous simulation with other apps
        
        ! Stress accumulation (daily)
        real(8) :: iqredwet_day = 0.0d0    ! T reduction due to oxygen stress
        real(8) :: iqreddry_day = 0.0d0    ! T reduction due to drought stress
        real(8) :: iqredsol_day = 0.0d0    ! T reduction due to salinity stress
        real(8) :: iqredfrs_day = 0.0d0    ! T reduction due to frost stress
        real(8) :: iptra_day = 0.0d0       ! Potential T since start of day
        
        ! Character variables
        character(len=11) :: date = ''     ! Current date string
        character(len=16) :: outfil = ''   ! Name of output file
        character(len=80) :: pathwork = '' ! Path to work directory
        character(len=80) :: project = ''  ! Name of project
        character(len=80) :: swpfile = ''  ! Name of main input file

        ! ------------------------------------------------------------------
        ! TimeControl persistent internals (formerly SAVE locals)
        ! Needed for multi-instance snapshot/restore while legacy code still
        ! calls TimeControl() with module variables.
        ! ------------------------------------------------------------------
        integer :: tc_datea(6) = 0
        integer :: tc_nextyear = 0
        integer :: tc_flprevious = 0
        logical :: tc_flTnext = .false.
        real(4) :: tc_fsec = 0.0
        real(8) :: tc_tchange = 0.0d0
        real(8) :: tc_dtEvent = 0.0d0
        real(8) :: tc_tEvent = 0.0d0
        real(8) :: tc_tcumold = 0.0d0
        real(8) :: tc_dtprevious = 0.0d0
        real(4) :: tc_tmptimestart = 0.0
        real(4) :: tc_tmptimeend = 0.0

        ! External/DLL exchange persistent state
        real(8) :: ex_tlast = 0.0d0
    end type time_state_t

    ! Soil state is provided by `soil_state_mod`.

    ! Atmosphere state is provided by `atmosphere_state_mod`.

    ! ===========================================================================
    ! Crop State
    ! ===========================================================================
    type :: crop_state_t
        ! Development stage
        real(8) :: dvs = 0.0d0                 ! Development stage (-)
        real(8) :: dvsend = 0.0d0              ! DVS at harvest
        real(8) :: tsum = 0.0d0                ! Temperature sum (°C)
        real(8) :: tsumgerm = 0.0d0            ! Temperature sum germination
        integer :: daycrop = 0                 ! Days since crop start
        integer :: daygrowth = 0               ! Days since grass management
        integer :: icrop = 0                   ! Current crop number
        
        ! Leaf area and canopy
        real(8) :: lai = 0.0d0                 ! Leaf area index
        real(8) :: laipot = 0.0d0              ! Potential LAI
        real(8) :: laimax = 0.0d0              ! Maximum LAI reached
        real(8) :: laiexp = 0.0d0              ! LAI in exponential stage
        real(8) :: gc = 0.0d0                  ! Ground cover (-)
        real(8) :: ch = 0.0d0                  ! Crop height (cm)
        real(8) :: cf = 0.0d0                  ! Crop factor (-)
        real(8) :: albedo = 0.0d0              ! Crop reflection
        
        ! Biomass partitions (actual)
        real(8) :: wlv = 0.0d0                 ! Leaf dry weight (kg/ha)
        real(8) :: wst = 0.0d0                 ! Stem dry weight
        real(8) :: wrt = 0.0d0                 ! Root dry weight
        real(8) :: wso = 0.0d0                 ! Storage organ dry weight
        real(8) :: cwdm = 0.0d0                ! Total dry weight
        real(8) :: tadw = 0.0d0                ! Above-ground dry weight
        
        ! Biomass partitions (potential)
        real(8) :: wlvpot = 0.0d0
        real(8) :: wstpot = 0.0d0
        real(8) :: wrtpot = 0.0d0
        real(8) :: wsopot = 0.0d0
        real(8) :: cwdmpot = 0.0d0
        real(8) :: tadwpot = 0.0d0
        
        ! Dead matter
        real(8) :: dwlv = 0.0d0                ! Dead leaves (kg/ha)
        real(8) :: dwst = 0.0d0                ! Dead stems
        real(8) :: dwrt = 0.0d0                ! Dead roots
        real(8) :: dwso = 0.0d0                ! Dead storage organs
        
        ! Assimilation and respiration
        real(8) :: gasst = 0.0d0               ! Gross assimilation
        real(8) :: gasstpot = 0.0d0            ! Potential gross assimilation
        real(8) :: mrest = 0.0d0               ! Maintenance respiration
        real(8) :: mrestpot = 0.0d0
        real(8) :: pgass = 0.0d0               ! Assimilation after N stress
        real(8) :: pgasspot = 0.0d0
        
        ! Rooting
        real(8) :: rd = 0.0d0                  ! Rooting depth (cm)
        real(8) :: rdpot = 0.0d0               ! Potential rooting depth
        real(8) :: rdi = 0.0d0                 ! Initial rooting depth
        real(8) :: rri = 0.0d0                 ! Max daily root increase
        real(8) :: rdc = 0.0d0                 ! Max crop rooting depth
        real(8) :: rdm = 0.0d0                 ! Max rooting depth
        real(8) :: rdmax = 0.0d0               ! Max profile rooting depth
        integer :: noddrz = 0                  ! Bottom node of root zone
        
        ! Root uptake parameters
        real(8) :: hlim1 = 0.0d0               ! Anaerobiosis point
        real(8) :: hlim2u = 0.0d0              ! Optimal point top
        real(8) :: hlim2l = 0.0d0              ! Optimal point bottom
        real(8) :: hlim3h = 0.0d0              ! Reduction start high Tpot
        real(8) :: hlim3l = 0.0d0              ! Reduction start low Tpot
        real(8) :: hlim4 = 0.0d0               ! Wilting point
        
        ! Stress factors
        real(8) :: reltr = 0.0d0               ! Relative transpiration factor
        real(8) :: alphacrit = 0.0d0           ! Critical stress index
        real(8) :: alpJvLier = 0.0d0           ! Jong van Lier reduction
        
        ! Grass specific
        real(8) :: tagp = 0.0d0                ! Total above-ground grass
        real(8) :: tagppot = 0.0d0
        real(8) :: tagpt = 0.0d0               ! Harvested grass
        real(8) :: mowrest = 0.0d0             ! Grass after mowing
        real(8) :: cuptgraz = 0.0d0            ! Cumulative grazing uptake
        integer :: iharvest = 0                ! Grass harvest number
        integer :: iseqgm = 0                  ! Grazing/mowing sequence counter
        
        ! Crop calendar
        real(8), allocatable :: cropstart(:)   ! Start dates per crop
        real(8), allocatable :: cropend(:)     ! End dates per crop
        
        ! Leaf age tracking
        real(8), allocatable :: lv(:)          ! Leaf weight by age
        real(8), allocatable :: lvage(:)       ! Leaf age
        real(8), allocatable :: sla(:)         ! Specific leaf area by age
        integer :: ilvold = 0                  ! Age of oldest leaf
        
        ! Flags
        logical :: flCropCalendar = .false.
        logical :: flCropEmergence = .false.
        logical :: flCropHarvest = .false.
        logical :: flanthesis = .false.
        logical :: flHarvest = .false.
        logical :: flHarvestDay = .false.
        logical :: flGrazing = .false.
        logical :: flCO2 = .false.
        
        ! CO2 correction factors
        real(8) :: fco2amax = 1.0d0
        real(8) :: fco2eff = 1.0d0
        real(8) :: fco2tra = 1.0d0
        
        ! Preparation/Sowing state
        logical :: flCropPrep = .false.
        logical :: flCropSow = .false.
        logical :: flCropGerm = .false.
        integer :: PrepDelay = 0
        integer :: SowDelay = 0
        integer :: DayGerm = 0
        
        ! Paths
        character(len=80) :: pathcrop = ''
        character(len=40), allocatable :: cropfil(:)
    end type crop_state_t

    ! ===========================================================================
    ! Irrigation State
    ! ===========================================================================
    type :: irrigation_state_t
        real(8) :: gird = 0.0d0                ! Gross irrigation depth
        real(8) :: nird = 0.0d0                ! Net irrigation depth
        real(8) :: cirrs = 0.0d0               ! Irrigation solute concentration
        real(8) :: cgird = 0.0d0               ! Cumulative gross irrigation
        real(8) :: cnird = 0.0d0               ! Cumulative net irrigation
        real(8) :: igird = 0.0d0               ! Intermediate gross irrigation
        real(8) :: inird = 0.0d0               ! Intermediate net irrigation
        
        integer :: irrigevent = 0              ! Current irrigation event
        integer :: nirri = 0                   ! Number of irrigation events
        integer :: isua = 0                    ! Irrigation type switch
        integer :: dayfix = 0                  ! Days since last irrigation
        
        logical :: flirrigate = .false.        ! Irrigation in simulation
        logical :: flheadirg = .false.         ! Print irrigation header
        
        ! SSDI (subsurface drip irrigation) state
        integer :: swssdi = 0                  ! Switch: SSDI active (0=no, 1=yes)
        integer, dimension(2) :: nod_ssdi = 0  ! Upper and lower nodes for SSDI
        integer :: ssdi_schedule = 0           ! Schedule type (0=fixed dates, 1=internal)
        integer :: ssdi_sched_type = 0         ! Internal schedule type (1=Tact/Tpot, 2=h, 3=theta)
        integer :: nod_ssdi_sensor = 0         ! Sensor node (if ssdi_sched_type > 1)
        real(8) :: ssdi_threshold = 0.0d0      ! Threshold value for scheduling
        real(8) :: ssdi_threshold_z = 0.0d0    ! Depth for threshold value
        real(8) :: ssdi_amount = 0.0d0         ! Amount of scheduled irrigation (cm)
        real(8) :: ssdi_appl_rate = 0.0d0      ! Application rate (cm/d)
        integer :: sw_interval = 0             ! Switch for minimum interval
        integer :: days_interval = 0           ! Minimum days between applications
        integer :: days_counter = 0            ! Days since previous application
        integer :: nirri_ssdi = 0              ! SSDI counter/entry point
        real(8), allocatable :: ssdi_date(:)   ! Fixed irrigation dates
        real(8), allocatable :: ssdi_rate_f(:) ! Fixed irrigation rates (cm/d)
        real(8), allocatable :: ssdi_amount_f(:) ! Fixed irrigation amounts (cm)
    end type irrigation_state_t

    ! ===========================================================================
    ! Tillage State
    ! ===========================================================================
    type :: tillage_state_t
        ! Switches and configuration
        integer :: swtill = 0                  ! Switch: 0=no tillage, 1=tillage
        integer :: i_n_model = 2               ! Switch for n-parameter treatment (1-3)
        integer :: iRedist = 2                 ! Redistribution type after MvG change
        
        ! Event counters
        integer :: Ntill = 0                   ! Number of tabulated tillage events
        integer :: iTill = 1                   ! Current tillage event index
        integer :: Ntypes = 0                  ! Number of tillage types
        integer :: MaxNumSoilHo = 0            ! Max soil horizons in tillage zone
        integer :: MaxNumSoilCP = 0            ! Max soil compartments in tillage zone
        
        ! Tillage depth
        real(8) :: Max_Z_tillage = 0.0d0       ! Max possible depth of tillage (cm)
        
        ! Tabulated input per tillage event (size: Ntill)
        real(8), allocatable :: Date_tillage(:)   ! Tillage dates
        real(8), allocatable :: Z_tillage(:)      ! Tillage depths (cm)
        real(8), allocatable :: I_tillage(:)      ! Tillage intensity (0-1)
        integer, allocatable :: Type_Tillage(:)   ! Tillage type index
        
        ! Tabulated input per tillage type (size: Ntypes)
        integer, allocatable :: iType_Tillage(:)  ! Tillage type identifier
        integer, allocatable :: iTT1(:)           ! First position in type table
        integer, allocatable :: iTT2(:)           ! Last position in type table
        real(8), allocatable :: TAB_Rho_tillage(:) ! Bulk density after tillage
        real(8), allocatable :: TAB_Rho_cons(:)    ! Consolidated bulk density
        real(8), allocatable :: TAB_K_R_cons(:)    ! Consolidation rate constant
        real(8), allocatable :: TAB_Rho_match(:)   ! Matching point density
        real(8), allocatable :: TAB_N_match(:)     ! Matching point n-value
        
        ! Per-layer state (size: NumLay)
        real(8), allocatable :: Rho_tillage(:)    ! Post-tillage bulk density per layer
        real(8), allocatable :: Rho_cons(:)       ! Consolidated density per layer
        real(8), allocatable :: Rho_last(:)       ! Previous density per layer
        real(8), allocatable :: K_R_cons(:)       ! Consolidation rate per layer
        real(8), allocatable :: Rho_match(:)      ! Matching point density per layer
        real(8), allocatable :: N_match(:)        ! Matching point n per layer
        real(8), allocatable :: Slope_match(:)    ! Slope at matching point per layer
        
        ! Water redistribution tracking
        real(8) :: sumDWC = 0.0d0              ! Sum of water content changes
        real(8) :: sumAvail1 = 0.0d0           ! Available pore space
        real(8) :: sumAvail2 = 0.0d0           ! Available water
    end type tillage_state_t

    ! ===========================================================================
    ! WOFOST Soil Nutrient State
    ! ===========================================================================
    type :: wofost_soil_state_t
        ! Interface variables (from wofost_soil_interface.f90)
        real(8) :: idwrt = 0.0d0               ! Increment dead root weight
        real(8) :: idwlv = 0.0d0               ! Increment dead leaf weight
        real(8) :: idwst = 0.0d0               ! Increment dead stem weight
        real(8) :: idwso = 0.0d0               ! Increment dead storage organ weight
        real(8) :: iNLOSSL = 0.0d0             ! N loss from leaves
        real(8) :: iNLOSSR = 0.0d0             ! N loss from roots
        real(8) :: iNLOSSS = 0.0d0             ! N loss from stems
        real(8) :: iNLOSSO = 0.0d0             ! N loss from storage organs
        real(8) :: NdemandSoil = 0.0d0         ! Total N demand (kg/ha N)
        real(8) :: NsupplySoil = 0.0d0         ! Total mineral N from soil (kg/ha N)
        real(8) :: Ndemand = 0.0d0             ! Total N demand (kg/m2 N)
        real(8) :: Nsupply = 0.0d0             ! Total mineral N supply (kg/m2 N)
        real(8) :: LaiCritNupt = 0.0d0         ! Critical LAI for N-uptake
        
        ! Soil water and temperature from SWAP
        real(8) :: dz_WSN = 0.0d0              ! Soil layer thickness
        real(8) :: WFrac_t = 0.0d0             ! Water fraction at time t
        real(8) :: WFrac_t0 = 0.0d0            ! Water fraction at time t0
        real(8) :: Wflux_out = 0.0d0           ! Water flux out
        real(8) :: Wflux_transp = 0.0d0        ! Transpiration flux
        real(8) :: Wflux_inBot = 0.0d0         ! Water flux in bottom
        real(8) :: Wflux_inTop = 0.0d0         ! Water flux in top
        real(8) :: Wflux_inLat = 0.0d0         ! Lateral water flux in
        real(8) :: SoilEvap = 0.0d0            ! Soil evaporation
        real(8) :: Temp_wsn = 0.0d0            ! Temperature for WSN
        real(8) :: t_WSNold = 0.0d0            ! Previous time for WSN
        
        ! Response function parameters
        real(8) :: Temp_ref = 0.0d0            ! Reference temperature
        real(8) :: WFrac_sat = 0.0d0           ! Saturated water fraction
        
        ! Organic matter model parameters (maxfn = 8)
        integer :: nf = 0                      ! Number of fractions
        real(8), dimension(8) :: RateconFOM_ref = 0.0d0  ! FOM rate constant ref
        real(8) :: RateconBio_ref = 0.0d0      ! Biomass rate constant ref
        real(8) :: RateconHum_ref = 0.0d0      ! Humus rate constant ref
        real(8) :: RateconHum_exp = 0.0d0      ! Humus rate constant exponent
        real(8) :: tstartHumexp = 0.0d0        ! Start time for Hum exponent
        real(8) :: t1900Soil = 0.0d0           ! Reference time 1900
        real(8), dimension(8) :: RateconFOM = 0.0d0      ! FOM rate constants
        real(8) :: RateconBio = 0.0d0          ! Biomass rate constant
        real(8) :: RateconHum = 0.0d0          ! Humus rate constant
        real(8), dimension(8) :: AsfaFOM_Bio = 0.0d0     ! Assimilation FOM->Bio
        real(8), dimension(8) :: AsfaFOM_Hum = 0.0d0     ! Assimilation FOM->Hum
        real(8) :: AsfaBio = 0.0d0             ! Assimilation Bio
        real(8) :: AsfaHum = 0.0d0             ! Assimilation Hum
        real(8), dimension(8) :: CFracFOM = 0.0d0        ! C fraction FOM
        real(8) :: CFracBio = 0.0d0            ! C fraction Bio
        real(8) :: CFracHum = 0.0d0            ! C fraction Hum
        real(8), dimension(8) :: NFracFOM = 0.0d0        ! N fraction FOM
        real(8) :: NFracBio = 0.0d0            ! N fraction Bio
        real(8) :: NFracHum = 0.0d0            ! N fraction Hum
        real(8) :: NFracFOMmin = 0.0d0         ! Min N fraction FOM
        real(8) :: NFracFOMmax = 0.0d0         ! Max N fraction FOM
        real(8) :: AsfaMin = 0.0d0             ! Min assimilation factor
        real(8) :: AsfaMax = 0.0d0             ! Max assimilation factor
        
        ! Organic matter state variables
        real(8), dimension(8) :: FOM_t0 = 0.0d0          ! FOM at t0
        real(8) :: Bio_t0 = 0.0d0              ! Biomass at t0
        real(8) :: Hum_t0 = 0.0d0              ! Humus at t0
        real(8), dimension(8) :: FOM_t = 0.0d0           ! FOM at t
        real(8) :: Bio_t = 0.0d0               ! Biomass at t
        real(8) :: Hum_t = 0.0d0               ! Humus at t
        real(8) :: Cdissi = 0.0d0              ! C dissimilation
        
        ! Mineral nitrogen model
        real(8) :: RateConNitrif_ref = 0.0d0   ! Nitrification rate ref
        real(8) :: RateConDenitr_ref = 0.0d0   ! Denitrification rate ref
        real(8) :: RateConNitrif = 0.0d0       ! Nitrification rate
        real(8) :: RateConDenitr = 0.0d0       ! Denitrification rate
        real(8) :: TCSF_N = 0.0d0              ! Temperature correction
        real(8) :: Nminer = 0.0d0              ! N mineralization
        real(8) :: Ratecon = 0.0d0             ! Rate constant
        
        ! Ammonium and nitrate concentrations
        real(8) :: DryBD = 0.0d0               ! Dry bulk density
        real(8) :: SorpCoef = 0.0d0            ! Sorption coefficient
        real(8) :: cNH4_t0 = 0.0d0             ! NH4 at t0
        real(8) :: cNH4_t = 0.0d0              ! NH4 at t
        real(8) :: cNH4_av = 0.0d0             ! NH4 average
        real(8) :: cNO3_t0 = 0.0d0             ! NO3 at t0
        real(8) :: cNO3_t = 0.0d0              ! NO3 at t
        real(8) :: cNO3_av = 0.0d0             ! NO3 average
        
        ! Boundary concentrations
        real(8) :: cNH4N_top = 0.0d0           ! NH4-N at top
        real(8) :: cNH4N_lat = 0.0d0           ! NH4-N lateral
        real(8) :: cNH4N_seep = 0.0d0          ! NH4-N seepage
        real(8) :: cNO3N_top = 0.0d0           ! NO3-N at top
        real(8) :: cNO3N_lat = 0.0d0           ! NO3-N lateral
        real(8) :: cNO3N_seep = 0.0d0          ! NO3-N seepage
        real(8) :: wsn_Cseep = 0.0d0           ! WOFOST seepage concentration
        real(8) :: wsn_Ctop = 0.0d0            ! WOFOST top concentration
        real(8) :: wsn_Clat = 0.0d0            ! WOFOST lateral concentration
        
        ! Time control
        real(8) :: dt_WSN = 0.0d0              ! Time step for WSN
        
        ! Response function parameters
        real(8) :: WFPSCrit = 0.0d0            ! Critical WFPS
        real(8) :: WFPScrit2 = 0.0d0           ! Second critical WFPS
        real(8) :: CdissiHalf = 0.0d0          ! Half C dissimilation
        real(8) :: WFPS = 0.0d0                ! Water-filled pore space
        real(8) :: red_T = 0.0d0               ! Temperature reduction
        real(8) :: red_W = 0.0d0               ! Water reduction
        real(8) :: red_W_Nit = 0.0d0           ! Water reduction nitrification
        real(8) :: red_W_Den = 0.0d0           ! Water reduction denitrification
        real(8) :: red_Resp = 0.0d0            ! Respiration reduction
        
        ! N supply tracking
        real(8) :: NsupplyNH4N = 0.0d0         ! NH4-N supply
        real(8) :: NsupplyNO3N = 0.0d0         ! NO3-N supply
        
        ! Old values for balance (scalars - summed across fractions)
        real(8) :: FOM_old = 0.0d0             ! Old FOM (summed)
        real(8) :: Bio_old = 0.0d0             ! Old biomass
        real(8) :: Hum_old = 0.0d0             ! Old humus
        real(8) :: NFOM_old = 0.0d0            ! Old N in FOM (summed)
        real(8) :: NBio_old = 0.0d0            ! Old N in biomass
        real(8) :: NHum_old = 0.0d0            ! Old N in humus
        real(8) :: NH4_old = 0.0d0             ! Old NH4
        real(8) :: NO3_old = 0.0d0             ! Old NO3
        
        ! Crop residue tracking
        real(8) :: iNLOSSL_1 = 0.0d0           ! N loss leaves (prev step)
        real(8) :: iNLOSSR_1 = 0.0d0           ! N loss roots (prev step)
        real(8) :: iNLOSSS_1 = 0.0d0           ! N loss stems (prev step)
        real(8) :: iNLOSSO_1 = 0.0d0           ! N loss storage organs (prev step)
        real(8) :: idwrt_1 = 0.0d0             ! Dead root (prev step)
        real(8) :: idwlv_1 = 0.0d0             ! Dead leaves (prev step)
        real(8) :: idwst_1 = 0.0d0             ! Dead stems (prev step)
        real(8) :: idwso_1 = 0.0d0             ! Dead storage organs (prev step)
        
        ! ANIMO crop_ext tracking
        real(8) :: Ntotuptake = 0.0d0          ! Total N uptake
        real(8) :: Ptotuptake = 0.0d0          ! Total P uptake
        real(8) :: DMcressur = 0.0d0           ! Crop residue surface DM
        real(8) :: Ncressurf = 0.0d0           ! Crop residue surface N
        real(8) :: Pcressurf = 0.0d0           ! Crop residue surface P
        real(8) :: DMcresbott = 0.0d0          ! Crop residue bottom DM
        real(8) :: Ncresbott = 0.0d0           ! Crop residue bottom N
        real(8) :: Pcresbott = 0.0d0           ! Crop residue bottom P
        
        ! Amendment state
        integer :: iAmendTime = 0              ! Current amendment time index
        integer :: isme = 0                    ! Soil management event index
        real(8) :: NH4N_volat = 0.0d0          ! NH4-N volatilization
        real(8) :: NH4N_amend = 0.0d0          ! NH4-N from amendment
        real(8) :: NO3N_amend = 0.0d0          ! NO3-N from amendment
        real(8) :: NH4N_cres = 0.0d0           ! NH4-N from crop residue
        real(8) :: NO3N_cres = 0.0d0           ! NO3-N from crop residue
        
        ! File unit numbers
        integer :: nut = -1                    ! Nutrient output file unit
        integer :: cropext = -1                ! Crop extension file unit
        
        ! Flags
        logical :: flCropExt = .false.         ! Crop extension file flag
    end type wofost_soil_state_t

    ! ===========================================================================
    ! Oxygen Stress State (from O2_pars module and OxygenStress subroutine)
    ! ===========================================================================
    type :: oxygenstress_state_t
        ! Root properties (from O2_pars module)
        real(8) :: w_root = 0.0d0              ! Dry weight per root length (kg/m)
        real(8) :: w_root_z0 = 0.0d0           ! Root weight at depth
        real(8) :: root_radius = 0.0d0         ! Root radius (m)
        
        ! Soil physical state
        real(8) :: soil_temp = 0.0d0           ! Soil temperature (K)
        real(8) :: sat_water_cont = 0.0d0      ! Saturated water content
        real(8) :: gas_filled_porosity = 0.0d0 ! Gas-filled porosity
        
        ! Diffusion coefficients
        real(8) :: d_o2inwater = 0.0d0         ! O2 diffusion in water
        real(8) :: d_root = 0.0d0              ! Diffusion in root
        real(8) :: d_soil = 0.0d0              ! Soil diffusion
        
        ! Soil properties
        real(8) :: perc_org_mat = 0.0d0        ! Organic matter percentage
        real(8) :: soil_density = 0.0d0        ! Soil density (kg/m3)
        real(8) :: depth = 0.0d0               ! Compartment thickness (m)
        
        ! Shape factors
        real(8) :: shape_factor_microbialr = 0.0d0
        
        ! Respiration state
        real(8) :: r_microbial_z0 = 0.0d0      ! Microbial respiration rate
        
        ! Water film / O2 transport
        real(8) :: waterfilm_thickness = 0.0d0
        real(8) :: bunsencoeff = 0.0d0
        
        ! O2 concentrations
        real(8) :: c_min_micro = 0.0d0         ! Min O2 for microbial resp
        real(8) :: c_macro = 0.0d0             ! Macropore O2 conc
        real(8) :: ctopnode = 0.0d0            ! Top node O2 conc
        
        ! Pre-calculated constants per node (from OxygenStress subroutine SAVE)
        real(8), allocatable :: d_soil_term1(:)
        real(8), allocatable :: d_soil_term2(:)
        real(8), allocatable :: gfp100(:)
        real(8), allocatable :: Capac_term(:)
        real(8), allocatable :: Nmin1(:)
        real(8), allocatable :: Mplus1(:)
        
        ! Initialization flags
        logical :: ini_stress = .true.         ! O2 stress initialization flag (from OxygenStress)
        logical :: initialized = .false.       ! Overall initialization flag
    end type oxygenstress_state_t

    ! Boundary state is provided by `boundary_state_mod`.

    ! Macropore state is provided by `macropore_state_mod`.

    ! ===========================================================================
    ! Snow State
    ! ===========================================================================
    type :: snow_state_t
        real(8) :: ssnow = 0.0d0               ! Snow water equivalent (cm)
        real(8) :: slw = 0.0d0                 ! Liquid water in snow
        real(8) :: gsnow = 0.0d0               ! Gross snow rate
        real(8) :: snrai = 0.0d0               ! Net rain on snow
        real(8) :: melt = 0.0d0                ! Melt rate
        real(8) :: subl = 0.0d0                ! Sublimation rate
        real(8) :: snowcoef = 0.0d0            ! Snow melt factor
        real(8) :: snowinco = 0.0d0            ! Snow at balance start
        
        ! Cumulative
        real(8) :: cgsnow = 0.0d0              ! Cumulative gross snow
        real(8) :: csnrai = 0.0d0              ! Cumulative net snow
        real(8) :: cmelt = 0.0d0               ! Cumulative melt
        real(8) :: csubl = 0.0d0               ! Cumulative sublimation
        
        ! Intermediate
        real(8) :: igsnow = 0.0d0
        real(8) :: isnrai = 0.0d0
        real(8) :: isubl = 0.0d0
        
        ! Temperature thresholds
        real(8) :: TePrRain = 0.0d0            ! All precip as rain above this
        real(8) :: TePrSnow = 0.0d0            ! All precip as snow below this
        
        ! Flags
        logical :: flsnow = .false.
    end type snow_state_t

    ! ===========================================================================
    ! I/O File Handles (Separate from simulation state)
    ! ===========================================================================
    type :: io_handles_t
        ! Input files
        integer :: swp_unit = -1               ! Main input file (.swp)
        integer :: met_unit = -1               ! Meteorology file
        integer :: crp_unit = -1               ! Crop file
        integer :: dra_unit = -1               ! Drainage file
        
        ! Output files
        integer :: log_unit = -1               ! Log file (.log)
        integer :: bal_unit = -1               ! Balance output (.bal)
        integer :: blc_unit = -1               ! Detailed balance (.blc)
        integer :: wba_unit = -1               ! Daily water balance (.wba)
        integer :: inc_unit = -1               ! Incremental balance (.inc)
        integer :: vap_unit = -1               ! Profile output (.vap)
        integer :: tem_unit = -1               ! Temperature output (.tem)
        integer :: str_unit = -1               ! Stress output (.str)
        integer :: crp_out_unit = -1           ! Crop output (.crp)
        integer :: irg_unit = -1               ! Irrigation output (.irg)
        integer :: snw_unit = -1               ! Snow output (.snw)
        integer :: sba_unit = -1               ! Solute balance (.sba)
        integer :: bma_unit = -1               ! Macropore balance (.bma)
        integer :: rot_unit = -1               ! Root extraction (.rot)
        integer :: drf_unit = -1               ! Drainage output
        integer :: swb_unit = -1               ! Surface water balance
        integer :: csv_unit = -1               ! CSV output
        
        ! External coupling files
        integer :: afo_unit = -1               ! Formatted for water quality
        integer :: aun_unit = -1               ! Unformatted for water quality
        
        ! State flags
        logical :: files_open = .false.
    end type io_handles_t

    ! ===========================================================================
    ! Top-Level State Container
    ! ===========================================================================
    type :: swap_state_t
        type(time_state_t)       :: time
        type(soil_state_t)       :: soil
        type(atmosphere_state_t) :: atm
        type(crop_state_t)       :: crop
        type(irrigation_state_t) :: irrig
        type(tillage_state_t)    :: tillage
        type(drainage_state_t)   :: drain
        type(boundary_state_t)   :: boundary
        type(macropore_state_t)  :: macro
        type(solute_state_t)     :: solute
        type(heat_state_t)       :: heat
        type(snow_state_t)       :: snow
        type(surfacewater_state_t) :: surfwater
        type(wofost_soil_state_t) :: wofost_soil
        type(oxygenstress_state_t) :: oxystress
        
        ! Model configuration (set once at initialization)
        integer :: numnod = 0                  ! Number of compartments
        integer :: numlay = 0                  ! Number of soil layers
        integer :: nrlevs = 0                  ! Number of drainage levels
        integer :: ncrop = 0                   ! Number of crops
        
        ! Initialization flag
        logical :: initialized = .false.
    end type swap_state_t

contains

    ! ===========================================================================
    ! Initialization Procedures
    ! ===========================================================================
    
    subroutine swap_state_init(state, numnod, numlay, nrlevs, ncrop)
        !> Initialize SWAP state with allocated arrays
        type(swap_state_t), intent(inout) :: state
        integer, intent(in) :: numnod          ! Number of compartments
        integer, intent(in) :: numlay          ! Number of soil layers
        integer, intent(in), optional :: nrlevs ! Number of drainage levels
        integer, intent(in), optional :: ncrop  ! Number of crops
        
        integer :: nlev, nc
        
        ! Set defaults for optional parameters
        nlev = 1
        nc = 1
        if (present(nrlevs)) nlev = nrlevs
        if (present(ncrop)) nc = ncrop
        
        ! Store configuration
        state%numnod = numnod
        state%numlay = numlay
        state%nrlevs = nlev
        state%ncrop = nc
        
        ! Initialize sub-states
        call soil_state_init(state%soil, numnod, numlay)

        call atmosphere_state_init(state%atm)
        
        call drainage_state_init(state%drain, nlev, numnod)
        
        call crop_state_init(state%crop, nc)
        
        call macropore_state_init(state%macro, numnod, MADM, MADR)
        
        call solute_state_init(state%solute, numnod, numlay, nlev)
        
        call heat_state_init(state%heat, numnod, numlay)
        
        call boundary_state_init(state%boundary)
        
        ! Initialize surfacewater state with reasonable defaults
        ! nmper=10 (management periods), mamte=100 (meteo entries), maowl=100 (water level entries)
        call surfacewater_state_init(state%surfwater, 10, 100, 100, nlev)
        
        call oxygenstress_state_init(state%oxystress, numnod)
        
        state%initialized = .true.
        
        call log_info('state_init', 'SWAP state initialized: ' // &
                      trim(to_str(numnod)) // ' nodes, ' // &
                      trim(to_str(numlay)) // ' layers')
        
    end subroutine swap_state_init
    
    subroutine crop_state_init(crop, ncrop)
        type(crop_state_t), intent(inout) :: crop
        integer, intent(in) :: ncrop
        
        allocate(crop%cropstart(ncrop))
        allocate(crop%cropend(ncrop))
        allocate(crop%cropfil(ncrop))
        allocate(crop%lv(366))
        allocate(crop%lvage(366))
        allocate(crop%sla(366))
        
        crop%cropstart = 0.0d0
        crop%cropend = 0.0d0
        crop%lv = 0.0d0
        crop%lvage = 0.0d0
        crop%sla = 0.0d0
        
    end subroutine crop_state_init
    
    subroutine io_handles_init(io)
        !> Initialize I/O handles to invalid state
        type(io_handles_t), intent(inout) :: io
        
        io%swp_unit = -1
        io%met_unit = -1
        io%crp_unit = -1
        io%dra_unit = -1
        io%log_unit = -1
        io%bal_unit = -1
        io%blc_unit = -1
        io%wba_unit = -1
        io%inc_unit = -1
        io%vap_unit = -1
        io%tem_unit = -1
        io%str_unit = -1
        io%crp_out_unit = -1
        io%irg_unit = -1
        io%snw_unit = -1
        io%sba_unit = -1
        io%bma_unit = -1
        io%rot_unit = -1
        io%drf_unit = -1
        io%swb_unit = -1
        io%csv_unit = -1
        io%afo_unit = -1
        io%aun_unit = -1
        io%files_open = .false.
        
    end subroutine io_handles_init
    
    !> @brief Finalize drainage state (deallocate arrays)
    subroutine drain_state_finalize(drain)
        type(drainage_state_t), intent(inout) :: drain

        call drainage_state_finalize(drain)

    end subroutine drain_state_finalize
    
    !> @brief Initialize irrigation state with allocated arrays
    subroutine irrigation_state_init(istate, maxirrig)
        type(irrigation_state_t), intent(inout) :: istate
        integer, intent(in) :: maxirrig          ! Maximum number of irrigation events
        
        ! Allocate SSDI arrays
        if (.not. allocated(istate%ssdi_date)) allocate(istate%ssdi_date(maxirrig))
        if (.not. allocated(istate%ssdi_rate_f)) allocate(istate%ssdi_rate_f(maxirrig))
        if (.not. allocated(istate%ssdi_amount_f)) allocate(istate%ssdi_amount_f(maxirrig))
        
        ! Initialize arrays to zero
        istate%ssdi_date = 0.0d0
        istate%ssdi_rate_f = 0.0d0
        istate%ssdi_amount_f = 0.0d0
        
    end subroutine irrigation_state_init
    
    !> @brief Finalize irrigation state (deallocate arrays)
    subroutine irrigation_state_finalize(istate)
        type(irrigation_state_t), intent(inout) :: istate
        
        if (allocated(istate%ssdi_date)) deallocate(istate%ssdi_date)
        if (allocated(istate%ssdi_rate_f)) deallocate(istate%ssdi_rate_f)
        if (allocated(istate%ssdi_amount_f)) deallocate(istate%ssdi_amount_f)
        
    end subroutine irrigation_state_finalize
    
    !> @brief Initialize tillage state with allocated arrays
    subroutine tillage_state_init(tstate, numlay, maxtill, maxtypes)
        type(tillage_state_t), intent(inout) :: tstate
        integer, intent(in) :: numlay             ! Number of soil layers
        integer, intent(in) :: maxtill            ! Maximum number of tillage events
        integer, intent(in) :: maxtypes           ! Maximum number of tillage types
        
        ! Allocate per-event arrays
        if (.not. allocated(tstate%Date_tillage)) allocate(tstate%Date_tillage(maxtill + 1))
        if (.not. allocated(tstate%Z_tillage)) allocate(tstate%Z_tillage(maxtill))
        if (.not. allocated(tstate%I_tillage)) allocate(tstate%I_tillage(maxtill))
        if (.not. allocated(tstate%Type_Tillage)) allocate(tstate%Type_Tillage(maxtill))
        
        ! Allocate per-type arrays
        if (.not. allocated(tstate%iType_Tillage)) allocate(tstate%iType_Tillage(maxtypes))
        if (.not. allocated(tstate%iTT1)) allocate(tstate%iTT1(maxtill))
        if (.not. allocated(tstate%iTT2)) allocate(tstate%iTT2(maxtill))
        if (.not. allocated(tstate%TAB_Rho_tillage)) allocate(tstate%TAB_Rho_tillage(maxtypes))
        if (.not. allocated(tstate%TAB_Rho_cons)) allocate(tstate%TAB_Rho_cons(maxtypes))
        if (.not. allocated(tstate%TAB_K_R_cons)) allocate(tstate%TAB_K_R_cons(maxtypes))
        if (.not. allocated(tstate%TAB_Rho_match)) allocate(tstate%TAB_Rho_match(maxtypes))
        if (.not. allocated(tstate%TAB_N_match)) allocate(tstate%TAB_N_match(maxtypes))
        
        ! Allocate per-layer arrays
        if (.not. allocated(tstate%Rho_tillage)) allocate(tstate%Rho_tillage(numlay))
        if (.not. allocated(tstate%Rho_cons)) allocate(tstate%Rho_cons(numlay))
        if (.not. allocated(tstate%Rho_last)) allocate(tstate%Rho_last(numlay))
        if (.not. allocated(tstate%K_R_cons)) allocate(tstate%K_R_cons(numlay))
        if (.not. allocated(tstate%Rho_match)) allocate(tstate%Rho_match(numlay))
        if (.not. allocated(tstate%N_match)) allocate(tstate%N_match(numlay))
        if (.not. allocated(tstate%Slope_match)) allocate(tstate%Slope_match(numlay))
        
        ! Initialize arrays to zero
        tstate%Date_tillage = 0.0d0
        tstate%Z_tillage = 0.0d0
        tstate%I_tillage = 0.0d0
        tstate%Type_Tillage = 0
        tstate%iType_Tillage = 0
        tstate%iTT1 = 0
        tstate%iTT2 = 0
        tstate%TAB_Rho_tillage = 0.0d0
        tstate%TAB_Rho_cons = 0.0d0
        tstate%TAB_K_R_cons = 0.0d0
        tstate%TAB_Rho_match = 0.0d0
        tstate%TAB_N_match = 0.0d0
        tstate%Rho_tillage = 0.0d0
        tstate%Rho_cons = 0.0d0
        tstate%Rho_last = 0.0d0
        tstate%K_R_cons = 0.0d0
        tstate%Rho_match = 0.0d0
        tstate%N_match = 0.0d0
        tstate%Slope_match = 0.0d0
        
    end subroutine tillage_state_init
    
    !> @brief Finalize tillage state (deallocate arrays)
    subroutine tillage_state_finalize(tstate)
        type(tillage_state_t), intent(inout) :: tstate
        
        ! Deallocate per-event arrays
        if (allocated(tstate%Date_tillage)) deallocate(tstate%Date_tillage)
        if (allocated(tstate%Z_tillage)) deallocate(tstate%Z_tillage)
        if (allocated(tstate%I_tillage)) deallocate(tstate%I_tillage)
        if (allocated(tstate%Type_Tillage)) deallocate(tstate%Type_Tillage)
        
        ! Deallocate per-type arrays
        if (allocated(tstate%iType_Tillage)) deallocate(tstate%iType_Tillage)
        if (allocated(tstate%iTT1)) deallocate(tstate%iTT1)
        if (allocated(tstate%iTT2)) deallocate(tstate%iTT2)
        if (allocated(tstate%TAB_Rho_tillage)) deallocate(tstate%TAB_Rho_tillage)
        if (allocated(tstate%TAB_Rho_cons)) deallocate(tstate%TAB_Rho_cons)
        if (allocated(tstate%TAB_K_R_cons)) deallocate(tstate%TAB_K_R_cons)
        if (allocated(tstate%TAB_Rho_match)) deallocate(tstate%TAB_Rho_match)
        if (allocated(tstate%TAB_N_match)) deallocate(tstate%TAB_N_match)
        
        ! Deallocate per-layer arrays
        if (allocated(tstate%Rho_tillage)) deallocate(tstate%Rho_tillage)
        if (allocated(tstate%Rho_cons)) deallocate(tstate%Rho_cons)
        if (allocated(tstate%Rho_last)) deallocate(tstate%Rho_last)
        if (allocated(tstate%K_R_cons)) deallocate(tstate%K_R_cons)
        if (allocated(tstate%Rho_match)) deallocate(tstate%Rho_match)
        if (allocated(tstate%N_match)) deallocate(tstate%N_match)
        if (allocated(tstate%Slope_match)) deallocate(tstate%Slope_match)
        
        ! Reset scalar values
        tstate%swtill = 0
        tstate%Ntill = 0
        tstate%iTill = 1
        tstate%Ntypes = 0
        
    end subroutine tillage_state_finalize
    
    !> @brief Initialize WOFOST soil state (no dynamic allocations needed)
    subroutine wofost_soil_state_init(wstate)
        type(wofost_soil_state_t), intent(inout) :: wstate
        
        ! Most fields have default initialization
        ! This subroutine is provided for completeness
        wstate%nut = -1
        wstate%cropext = -1
        wstate%flCropExt = .false.
        
    end subroutine wofost_soil_state_init
    
    !> @brief Finalize WOFOST soil state (no dynamic allocations)
    subroutine wofost_soil_state_finalize(wstate)
        type(wofost_soil_state_t), intent(inout) :: wstate
        
        ! No dynamic allocations to deallocate
        ! Reset file units
        wstate%nut = -1
        wstate%cropext = -1
        
    end subroutine wofost_soil_state_finalize
    
    !> @brief Initialize oxygen stress state
    subroutine oxygenstress_state_init(ostate, n_nod)
        type(oxygenstress_state_t), intent(inout) :: ostate
        integer, intent(in) :: n_nod
        
        ! Allocate per-node arrays
        if (.not. allocated(ostate%d_soil_term1)) allocate(ostate%d_soil_term1(n_nod))
        if (.not. allocated(ostate%d_soil_term2)) allocate(ostate%d_soil_term2(n_nod))
        if (.not. allocated(ostate%gfp100)) allocate(ostate%gfp100(n_nod))
        if (.not. allocated(ostate%Capac_term)) allocate(ostate%Capac_term(n_nod))
        if (.not. allocated(ostate%Nmin1)) allocate(ostate%Nmin1(n_nod))
        if (.not. allocated(ostate%Mplus1)) allocate(ostate%Mplus1(n_nod))
        
        ! Initialize arrays
        ostate%d_soil_term1 = 0.0d0
        ostate%d_soil_term2 = 0.0d0
        ostate%gfp100 = 0.0d0
        ostate%Capac_term = 0.0d0
        ostate%Nmin1 = 0.0d0
        ostate%Mplus1 = 0.0d0
        
        ostate%initialized = .false.
        
    end subroutine oxygenstress_state_init
    
    !> @brief Finalize oxygen stress state
    subroutine oxygenstress_state_finalize(ostate)
        type(oxygenstress_state_t), intent(inout) :: ostate
        
        if (allocated(ostate%d_soil_term1)) deallocate(ostate%d_soil_term1)
        if (allocated(ostate%d_soil_term2)) deallocate(ostate%d_soil_term2)
        if (allocated(ostate%gfp100)) deallocate(ostate%gfp100)
        if (allocated(ostate%Capac_term)) deallocate(ostate%Capac_term)
        if (allocated(ostate%Nmin1)) deallocate(ostate%Nmin1)
        if (allocated(ostate%Mplus1)) deallocate(ostate%Mplus1)
        
        ostate%initialized = .false.
        
    end subroutine oxygenstress_state_finalize
    
    !> @brief Finalize (deallocate) SWAP state
    !> @param state The state to finalize
    subroutine swap_state_finalize(state)
        type(swap_state_t), intent(inout) :: state
        
        call soil_state_finalize(state%soil)
        call atmosphere_state_finalize(state%atm)
        call drain_state_finalize(state%drain)
        call surfacewater_state_finalize(state%surfwater)
        call boundary_state_finalize(state%boundary)
        call macropore_state_finalize(state%macro)
        call solute_state_finalize(state%solute)
        call heat_state_finalize(state%heat)
        call irrigation_state_finalize(state%irrig)
        call tillage_state_finalize(state%tillage)
        call wofost_soil_state_finalize(state%wofost_soil)
        call oxygenstress_state_finalize(state%oxystress)
        
        state%initialized = .false.
        
    end subroutine swap_state_finalize

end module swap_state_mod
