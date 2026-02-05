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
    use swap_log, only: log_info, log_error, to_str
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
    public :: irrigation_state_init
    public :: tillage_state_init
    public :: wofost_soil_state_init
    public :: oxygenstress_state_init
    
    ! Finalization procedures
    public :: swap_state_finalize
    public :: drain_state_finalize
    public :: surfacewater_state_finalize
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
    end type time_state_t

    ! ===========================================================================
    ! Soil Water State
    ! ===========================================================================
    type :: soil_state_t
        ! Primary state variables (node-based)
        real(8), allocatable :: h(:)           ! Pressure head (cm)
        real(8), allocatable :: hm1(:)         ! Pressure head at previous time
        real(8), allocatable :: theta(:)       ! Volumetric water content (-)
        real(8), allocatable :: thetm1(:)      ! Water content at previous time
        real(8), allocatable :: k(:)           ! Hydraulic conductivity (cm/d)
        real(8), allocatable :: kmean(:)       ! Mean K at interfaces
        real(8), allocatable :: dimoca(:)      ! Differential moisture capacity
        
        ! Geometry
        real(8), allocatable :: z(:)           ! Node depths (cm)
        real(8), allocatable :: dz(:)          ! Compartment thickness (cm)
        real(8), allocatable :: ztopcp(:)      ! Top of compartment
        real(8), allocatable :: zbotcp(:)      ! Bottom of compartment
        real(8), allocatable :: disnod(:)      ! Distance between nodes
        
        ! Soil properties per node
        real(8), allocatable :: thetar(:)      ! Residual water content
        real(8), allocatable :: thetas(:)      ! Saturated water content
        real(8), allocatable :: hroot(:)       ! Pressure head at root-soil interface
        real(8), allocatable :: twilt(:)       ! Wilting point pressure head
        real(8), allocatable :: rfcp(:)        ! Frost reduction factor
        
        ! Soil layer properties
        real(8), allocatable :: bdens(:)       ! Bulk density (g/cm3)
        real(8), allocatable :: cofani(:)      ! Anisotropy coefficient
        real(8), allocatable :: ksatfit(:)     ! Saturated K fitted
        real(8), allocatable :: ksatexm(:)     ! Saturated K measured
        real(8), allocatable :: thetsl(:)      ! Saturated water content per layer
        
        ! Van Genuchten parameters
        real(8), allocatable :: paramvg(:,:)   ! VG parameters (21,maho)
        real(8), allocatable :: cofgen(:,:)    ! Adjusted VG parameters (21,macp)
        
        ! Fluxes
        real(8), allocatable :: q(:)           ! Water flux between compartments
        real(8), allocatable :: evp(:)         ! Evaporation flux
        real(8), allocatable :: qrot(:)        ! Root extraction flux
        real(8), allocatable :: qpotrot(:)     ! Potential root extraction
        real(8), allocatable :: qssdi(:)       ! Subsurface drip irrigation
        
        ! Stress reduction fluxes
        real(8), allocatable :: qredwet(:)     ! Reduction due to wet conditions
        real(8), allocatable :: qreddry(:)     ! Reduction due to dry conditions
        real(8), allocatable :: qredsol(:)     ! Reduction due to salinity
        real(8), allocatable :: qredfrs(:)     ! Reduction due to frost
        
        ! Intermediate storage
        real(8), allocatable :: inq(:)         ! Intermediate water flux
        real(8), allocatable :: inqrot(:)      ! Intermediate root extraction
        real(8), allocatable :: inqssdi(:)     ! Intermediate SSDI
        real(8), allocatable :: ithetabeg(:)   ! Water content at period start
        
        ! Groundwater
        real(8) :: gwl = 0.0d0                 ! Groundwater level (cm)
        real(8) :: gwlm1 = 0.0d0               ! GWL at previous time
        real(8) :: gwli = 0.0d0                ! Initial GWL
        real(8) :: gwlinp = 0.0d0              ! Prescribed GWL
        real(8) :: pegwl = 0.0d0               ! Perched GWL
        real(8) :: deepgw = 0.0d0              ! Hydraulic head in aquifer
        integer :: nodgwl = 0                  ! Node above GWL
        integer :: npegwl = 0                  ! Node above perched GWL
        integer :: bpegwl = 0                  ! Node at bottom perched GWL
        
        ! Surface water/ponding
        real(8) :: pond = 0.0d0                ! Ponding depth (cm)
        real(8) :: pondm1 = 0.0d0              ! Ponding at previous time
        real(8) :: pondmx = 0.0d0              ! Maximum ponding
        real(8) :: pondini = 0.0d0             ! Initial ponding
        real(8) :: hsurf = 0.0d0               ! Pressure head at surface
        
        ! Bottom boundary
        real(8) :: qbot = 0.0d0                ! Bottom flux
        real(8) :: hbot = 0.0d0                ! Bottom pressure head
        
        ! Top boundary
        real(8) :: qtop = 0.0d0                ! Top flux
        logical :: ftoph = .false.             ! Prescribed pressure head at top
        
        ! Cumulative fluxes
        real(8) :: cqbot = 0.0d0               ! Cumulative bottom flux
        real(8) :: cqbotdo = 0.0d0             ! Cumulative downward bottom
        real(8) :: cqbotup = 0.0d0             ! Cumulative upward bottom
        real(8) :: cqtdo = 0.0d0               ! Cumulative downward top
        real(8) :: cqtup = 0.0d0               ! Cumulative upward top
        real(8) :: cqrot = 0.0d0               ! Cumulative root extraction
        real(8) :: cqdra = 0.0d0               ! Cumulative drainage
        real(8) :: crunoff = 0.0d0             ! Cumulative runoff
        real(8) :: crunon = 0.0d0              ! Cumulative runon
        
        ! Intermediate fluxes
        real(8) :: iqbot = 0.0d0
        real(8) :: iqrot = 0.0d0
        real(8) :: iqdra = 0.0d0
        real(8) :: iruno = 0.0d0
        real(8) :: irunon = 0.0d0
        real(8) :: iqssdi = 0.0d0
        real(8) :: ipondbeg = 0.0d0
        
        ! Storage
        real(8) :: volact = 0.0d0              ! Current water storage
        real(8) :: volini = 0.0d0              ! Initial water storage
        real(8) :: volm1 = 0.0d0               ! Previous water storage
        real(8) :: wbalance = 0.0d0            ! Water balance error
        
        ! Iteration control
        integer :: numbit = 0                  ! Iteration number
        integer :: msteps = 0                  ! Max iteration steps per day
        real(8) :: CritDevh1Cp = 0.0d0         ! Convergence criterion h relative
        real(8) :: CritDevh2Cp = 0.0d0         ! Convergence criterion h absolute
        real(8) :: CritDevMasBal = 0.0d0       ! Max water balance error
        real(8) :: gwlconv = 0.0d0             ! GWL convergence criterion
        
        ! Discretization info
        integer :: numnod = 0                  ! Number of nodes
        integer :: numlay = 0                  ! Number of soil layers
        integer :: nsublay = 0                 ! Number of sublayers
        integer, allocatable :: layer(:)       ! Layer number per compartment
        integer, allocatable :: botcom(:)      ! Bottom compartment per layer
        integer, allocatable :: nod1lay(:)     ! First node per layer
        
        ! Hysteresis
        integer :: swhyst = 0                  ! Hysteresis switch
        integer, allocatable :: indeks(:)      ! Wetting/drying curve index
        real(8) :: tau = 0.0d0                 ! Min h difference for transition
        
        ! Soil evaporation reduction
        real(8) :: saev = 0.0d0                ! Cumulative actual E (Boesten)
        real(8) :: spev = 0.0d0                ! Cumulative potential E (Boesten)
        real(8) :: ldwet = 0.0d0               ! Dry period length (Black)
        real(8) :: cofred = 0.0d0              ! Reduction coefficient
        
        ! Headcalc iteration tracking (from headcalc.f90 SAVE variables)
        logical :: flwarn_hc = .true.          ! Warning flag for headcalc
        integer :: iwarn_hc = 0                ! Warning counter
        integer :: nstep_hc = 0                ! Step counter in headcalc
        
        ! Flags
        logical :: FlRunoff = .false.
        logical :: fldrain = .false.
        logical :: flrunon = .false.
        logical :: fllowgwl = .false.
    end type soil_state_t

    ! ===========================================================================
    ! Atmosphere/Meteorology State
    ! ===========================================================================
    type :: atmosphere_state_t
        ! Current meteorological values
        real(8) :: tav = 0.0d0                 ! Average air temperature (°C)
        real(8) :: tavd = 0.0d0                ! Average daytime temperature
        real(8) :: tmn = 0.0d0                 ! Minimum temperature
        real(8) :: tmx = 0.0d0                 ! Maximum temperature
        real(8) :: tmnr = 0.0d0                ! 7-day average Tmin
        real(8) :: rad = 0.0d0                 ! Global solar radiation (J/m2/d)
        real(8) :: rh = 0.0d0                  ! Relative humidity (-)
        real(8) :: lat = 0.0d0                 ! Latitude (degrees)
        real(8) :: alt = 0.0d0                 ! Altitude (m)
        real(8) :: daylp = 0.0d0               ! Photoperiodic daylength (hours)
        
        ! Precipitation
        real(8) :: grai = 0.0d0                ! Gross rain flux (cm/d)
        real(8) :: graidt = 0.0d0              ! Gross precip during timestep
        real(8) :: nraida = 0.0d0              ! Net rain daily average
        real(8) :: nraidt = 0.0d0              ! Net rain during timestep
        real(8) :: finterception = 0.0d0       ! Net/gross rain ratio
        
        ! Evapotranspiration
        real(8) :: peva = 0.0d0                ! Potential soil evaporation (cm/d)
        real(8) :: pevaday = 0.0d0             ! Potential E of one day
        real(8) :: ptra = 0.0d0                ! Potential transpiration (cm/d)
        real(8) :: ptraday = 0.0d0             ! Potential T of one day
        real(8) :: tra = 0.0d0                 ! Actual transpiration (cm/d)
        real(8) :: reva = 0.0d0                ! Actual soil evaporation
        real(8) :: atmdem = 0.0d0              ! Atmospheric demand
        real(8) :: es0 = 0.0d0                 ! Potential E wet bare soil
        real(8) :: et0 = 0.0d0                 ! Potential T dry crop
        real(8) :: ew0 = 0.0d0                 ! Potential T wet crop
        
        ! Interception
        real(8) :: aintcdt = 0.0d0             ! Interception during timestep
        real(8) :: sicact = 0.0d0              ! Water stored on canopy
        real(8) :: siccapact = 0.0d0           ! Interception capacity
        
        ! Cumulative values
        real(8) :: cgrai = 0.0d0               ! Cumulative gross precipitation
        real(8) :: cnrai = 0.0d0               ! Cumulative net precipitation
        real(8) :: caintc = 0.0d0              ! Cumulative interception
        real(8) :: cevap = 0.0d0               ! Cumulative actual evaporation
        real(8) :: cpeva = 0.0d0               ! Cumulative potential evaporation
        real(8) :: cptra = 0.0d0               ! Cumulative potential transpiration
        
        ! Intermediate values
        real(8) :: inrai = 0.0d0
        real(8) :: igrai = 0.0d0
        real(8) :: iintc = 0.0d0
        real(8) :: ievap = 0.0d0
        real(8) :: ipeva = 0.0d0
        real(8) :: iptra = 0.0d0
        real(8) :: ies0 = 0.0d0
        real(8) :: iet0 = 0.0d0
        real(8) :: iew0 = 0.0d0
        
        ! Meteo file reading state
        integer :: daymeteo = 0
        integer :: yearmeteo = 0
        integer :: daynrfirst = 0
        integer :: daynrlast = 0
        integer :: wrecord = 0
        integer :: rainrec = 0
        real(8) :: timjan1 = 0.0d0
        real(8) :: metperiod = 0.0d0
        real(8) :: dtEventRain = 0.0d0
        
        ! Minimum temperature history (for 7-day average)
        real(8) :: atmin7(7) = 0.0d0
        
        ! Flags
        logical :: fletsine = .false.
        logical :: flmeteodt = .false.
        logical :: flmetdetail = .false.
        logical :: flrainintens = .false.
        logical :: flupdmetdet = .false.
        
        ! ETSine sub-daily state (from meteodt.f90 SAVE)
        real(8) :: tsunrise = 0.0d0            ! Time of sunrise (fraction of day)
        real(8) :: tsunset = 0.0d0             ! Time of sunset (fraction of day)
        
        ! CN runoff method state (from meteoday.f90 SAVE)
        integer :: nod10_cn = 0                ! Node at -10cm for CN method
        integer :: icn_atm = 0                 ! Current position in CN time table
        real(8) :: z10_cn = 0.0d0              ! Depth to node 10 for CN method
        
        ! Additional evaporation parameters
        real(8) :: empreva = 0.0d0             ! Reduced soil evaporation flux (L/T)
        real(8) :: fprecnosnow = 0.0d0         ! Ratio rain (excl snow) / gross rain
        
        ! Paths and files
        character(len=200) :: metfil = ''
        character(len=200) :: rainfil = ''
        character(len=80) :: pathatm = ''
    end type atmosphere_state_t

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

    ! ===========================================================================
    ! Drainage State
    ! ===========================================================================
    type :: drainage_state_t
        ! Drainage fluxes per level
        real(8), allocatable :: qdrain(:)       ! Drainage flux per level (cm/d)
        real(8), allocatable :: cqdrain(:)      ! Cumulative drainage per level
        real(8), allocatable :: cqdrainin(:)    ! Cumulative infiltration per level
        real(8), allocatable :: cqdrainout(:)   ! Cumulative drainage out per level
        real(8), allocatable :: drainl(:)       ! Drainage level depth
        
        ! Spatial drainage arrays
        real(8), allocatable :: qdra(:,:)       ! Drainage flux (level,node)
        real(8), allocatable :: inqdra(:,:)     ! Intermediate drainage (level,node)
        real(8), allocatable :: inqdra_in(:,:)  ! Intermediate infiltration flux per level/node
        real(8), allocatable :: inqdra_out(:,:) ! Intermediate drainage out per level/node
        real(8), allocatable :: qdraincomp(:)   ! Drainage per compartment
        
        ! Totals
        real(8) :: qdrtot = 0.0d0              ! Total drainage flux
        real(8) :: iqdra = 0.0d0               ! Intermediate drainage
        real(8) :: cqdra = 0.0d0               ! Cumulative total lateral drainage
        
        ! Parameters
        integer :: nrlevs = 0                  ! Number of drainage levels
        integer :: nrpri = 0                   ! Number of primary drainage levels
        integer :: dramet = 0                  ! Drainage method switch
        integer :: swdivd = 0                  ! Distribution of drainage over profile
        integer :: swdislay = 0                ! Discharge layer option
        real(8) :: basegw = 0.0d0              ! Impervious layer depth
        real(8) :: entres = 0.0d0              ! Drain entry resistance
        real(8) :: shape = 0.0d0               ! Shape factor
        
        ! Resistances and geometry per level
        real(8), allocatable :: drares(:)      ! Drainage resistance
        real(8), allocatable :: infres(:)      ! Infiltration resistance
        real(8), allocatable :: L(:)           ! Drain spacing
        real(8), allocatable :: wetper(:)      ! Wet perimeter
        real(8), allocatable :: zbotdr(:)      ! Drain bottom depth
        real(8), allocatable :: rdrain(:)      ! Drainage resistance per level
        real(8), allocatable :: rinfi(:)       ! Infiltration resistance per level
        real(8), allocatable :: rentry(:)      ! Entry resistance per level
        real(8), allocatable :: rexit(:)       ! Exit resistance per level
        real(8), allocatable :: gwlinf(:)      ! GWL below which no infiltration
        real(8), allocatable :: widthr(:)      ! Width of drain/channel
        real(8), allocatable :: taludr(:)      ! Talus slope of drain
        integer, allocatable :: swallo(:)      ! Allow drainage/infiltration per level
        integer, allocatable :: swdtyp(:)      ! Drainage type per level (0=channel, 1=tube, 2=interflow)
        integer, allocatable :: swtopdislay(:) ! Top of discharge layer option
        real(8), allocatable :: zTopDisLay(:)  ! Top of discharge layer
        real(8), allocatable :: fTopDisLay(:)  ! Fraction for discharge layer top
        
        ! Interflow parameters
        real(8) :: cofintfl = 0.0d0            ! Interflow coefficient
        real(8) :: expintfl = 0.0d0            ! Interflow exponent
        integer :: swnrsrf = 0                 ! Interflow switch
        integer :: SwTopnrsrf = 0              ! Interflow top switch
        real(8) :: rsurfdeep = 0.0d0           ! Deep interflow resistance
        real(8) :: rsurfshallow = 0.0d0        ! Shallow interflow resistance
        real(8) :: FacDpthInf = 0.0d0          ! Factor for infiltration depth
        integer :: Swdivdinf = 0               ! Infiltration distribution switch
        
        ! Flags
        logical :: fldrain = .false.
        
        ! Paths
        character(len=16) :: drfil = ''
        character(len=80) :: pathdrain = ''
    end type drainage_state_t

    ! ===========================================================================
    ! Boundary Conditions State
    ! ===========================================================================
    type :: boundary_state_t
        ! Bottom boundary - configuration
        integer :: swbotb = 0                  ! Bottom BC type (1-8)
        integer :: swbotb3Impl = 0             ! Implicit solution switch for BC type 3
        integer :: SwBotb3ResVert = 0          ! Suppress vertical resistance for BC type 3
        integer :: swqhbot = 0                 ! Flux-GWL relationship type
        integer :: swcofqhc = 0                ! Additional flux switch
        integer :: sw2 = 0                     ! Sub-switch for BC type 2
        integer :: sw3 = 0                     ! Sub-switch for BC type 3
        integer :: sw4 = 0                     ! Sub-switch for BC type 4 (extra flux)
        
        ! Bottom boundary - values
        real(8) :: qbot = 0.0d0                ! Bottom flux
        real(8) :: qbot_nonfrozen = 0.0d0      ! Bottom flux for non-frozen soil
        real(8) :: hbot = 0.0d0                ! Bottom pressure head
        real(8) :: iqbot = 0.0d0               ! Intermediate bottom flux
        real(8) :: cqbot = 0.0d0               ! Cumulative bottom flux
        real(8) :: cqbotdo = 0.0d0             ! Cumulative downward bottom flux
        real(8) :: cqbotup = 0.0d0             ! Cumulative upward bottom flux
        real(8) :: deepgw = 0.0d0              ! Hydraulic head in aquifer
        
        ! Aquifer parameters
        real(8) :: aqave = 0.0d0               ! Average aquifer head
        real(8) :: aqamp = 0.0d0               ! Aquifer head amplitude
        real(8) :: aqper = 0.0d0               ! Period of prescribed sine wave
        real(8) :: aqtmax = 0.0d0              ! Time with maximum head
        real(8) :: rimlay = 0.0d0              ! Aquitard resistance
        real(8) :: hdrain = 0.0d0              ! Mean drainage level
        real(8) :: shape = 0.0d0               ! Shape factor
        
        ! Sine function parameters for bottom flux
        real(8) :: sinave = 0.0d0              ! Average bottom flux
        real(8) :: sinamp = 0.0d0              ! Amplitude of bottom flux
        real(8) :: sinmax = 0.0d0              ! Time of maximum flux
        
        ! Flux-head relationships
        real(8) :: cofqha = 0.0d0              ! Coefficient A (exp function)
        real(8) :: cofqhb = 0.0d0              ! Coefficient B
        real(8) :: cofqhc = 0.0d0              ! Coefficient C
        
        ! Lysimeter parameters
        real(8) :: hplate = 0.0d0              ! Pressure head of ceramic plate
        
        ! Prescribed tables
        real(8), allocatable :: gwltab(:)      ! GWL vs time
        real(8), allocatable :: haqtab(:)      ! Aquifer head vs time
        real(8), allocatable :: qbotab(:)      ! Bottom flux vs time
        real(8), allocatable :: hbotab(:)      ! Bottom h vs time
        
        ! Top boundary - configuration
        integer :: swpondmx = 0                ! Time-dependent pondmx switch
        integer :: swredu = 0                  ! Soil evaporation reduction switch
        
        ! Top boundary - values
        real(8) :: pondmx = 0.0d0              ! Max ponding depth
        real(8) :: hatm = 0.0d0                ! Atmospheric pressure head
        real(8) :: hsurf = 0.0d0               ! Surface pressure head
        real(8) :: rsro = 0.0d0                ! Runoff resistance
        real(8) :: rsroexp = 0.0d0             ! Runoff exponent
        real(8) :: runon = 0.0d0               ! Runon flux
        real(8) :: runots = 0.0d0              ! Runoff during timestep
        real(8) :: crunoff = 0.0d0             ! Cumulative runoff
        real(8) :: crunon = 0.0d0              ! Cumulative runon
        real(8) :: iruno = 0.0d0               ! Intermediate runoff
        real(8) :: irunon = 0.0d0              ! Intermediate runon
        
        ! Top boundary - surface/ponding
        real(8) :: qtop = 0.0d0                ! Surface flux
        real(8) :: q0 = 0.0d0                  ! Net surface flux (precip-evap)
        real(8) :: h0max = 0.0d0               ! Max ponding without runoff
        real(8) :: k1max = 0.0d0               ! Max conductivity at surface
        real(8) :: QMpLatSs = 0.0d0            ! Lateral inflow to macropores at surface
        
        ! Runon table
        real(8), allocatable :: runonarr(:)    ! Daily runon values
        real(8), allocatable :: pondmxtab(:)   ! Time-dependent pondmx table
        
        ! Inundation
        real(8) :: cinund = 0.0d0              ! Cumulative inundation
        
        ! Flags
        logical :: FlRunoff = .false.          ! Runoff potential possible
        logical :: flrunon = .false.           ! Runon exists
        logical :: ftoph = .false.             ! Pressure head prescribed at surface
    end type boundary_state_t

    ! ===========================================================================
    ! Macropore State
    ! ===========================================================================
    type :: macropore_state_t
        ! Domain water storage
        real(8) :: VlMp = 0.0d0                ! Total macropore volume
        real(8) :: VlMpDm1 = 0.0d0             ! MB domain volume
        real(8) :: VlMpDm2 = 0.0d0             ! IC domain volume
        real(8) :: WaSrDm1 = 0.0d0             ! MB water storage
        real(8) :: WaSrDm2 = 0.0d0             ! IC water storage
        real(8) :: WaSrDm1Ini = 0.0d0          ! Initial MB storage
        real(8) :: WaSrDm2Ini = 0.0d0          ! Initial IC storage
        real(8) :: WaLevDm1 = 0.0d0            ! Water level in MB
        
        ! Fluxes
        real(8) :: QMaPo = 0.0d0               ! Matrix-macropore exchange
        real(8) :: QRapDra = 0.0d0             ! Rapid drainage flux
        real(8) :: QMpLatSs = 0.0d0            ! Lateral inflow at surface
        real(8) :: QInTopLatDm1 = 0.0d0        ! Top lateral inflow MB
        real(8) :: QInTopLatDm2 = 0.0d0        ! Top lateral inflow IC
        real(8) :: QInTopVrtDm1 = 0.0d0        ! Top vertical inflow MB
        real(8) :: QInTopVrtDm2 = 0.0d0        ! Top vertical inflow IC
        
        ! Cumulative fluxes
        real(8) :: cQMpLatSs = 0.0d0
        real(8) :: cQMpOutDrRap = 0.0d0
        real(8) :: cQMpInMtxSatDm1 = 0.0d0
        real(8) :: cQMpInMtxSatDm2 = 0.0d0
        real(8) :: cQMpOutMtxUnsDm1 = 0.0d0
        real(8) :: cQMpOutMtxUnsDm2 = 0.0d0
        real(8) :: cQMpInIntSatDm1 = 0.0d0
        real(8) :: cQMpInIntSatDm2 = 0.0d0
        real(8) :: cQMpInTopLatDm1 = 0.0d0
        real(8) :: cQMpInTopLatDm2 = 0.0d0
        real(8) :: cQMpInTopVrtDm1 = 0.0d0
        real(8) :: cQMpInTopVrtDm2 = 0.0d0
        real(8) :: cQMpOutMtxSatDm1 = 0.0d0
        real(8) :: cQMpOutMtxSatDm2 = 0.0d0
        
        ! Incremental fluxes
        real(8) :: iQMpOutDrRap = 0.0d0
        real(8) :: iQInTopLatDm1 = 0.0d0
        real(8) :: iQInTopLatDm2 = 0.0d0
        real(8) :: iQInTopVrtDm1 = 0.0d0
        real(8) :: iQInTopVrtDm2 = 0.0d0
        real(8) :: IWaSrDm1Beg = 0.0d0
        real(8) :: IWaSrDm2Beg = 0.0d0
        
        ! Per-compartment arrays
        real(8), allocatable :: DiPoCp(:)      ! Polygon diameter
        real(8), allocatable :: FrArMtrx(:)    ! Matrix area fraction
        real(8), allocatable :: VlMpDyCp(:)    ! Dynamic macropore volume
        real(8), allocatable :: VlMpStCp(:)    ! Static macropore volume
        real(8), allocatable :: VlMpStDm1(:)   ! Static volume domain 1
        real(8), allocatable :: VlMpStDm2(:)   ! Static volume domain 2
        real(8), allocatable :: QExcMpMtx(:)   ! Exchange flux
        real(8), allocatable :: SubsidCp(:)    ! Vertical subsidence
        real(8), allocatable :: PpDmCp(:,:)    ! Domain proportion
        real(8), allocatable :: dFdhMp(:)      ! Contribution to derivative
        real(8), allocatable :: iQOutDrRapCp(:)    ! Incremental rapid drainage per comp
        real(8), allocatable :: iQExcMtxDm1Cp(:)   ! Incremental exchange domain 1
        real(8), allocatable :: iQExcMtxDm2Cp(:)   ! Incremental exchange domain 2
        real(8), allocatable :: IAvFrMpWlWtDm1(:)  ! Avg wet wall fraction MB
        real(8), allocatable :: IAvFrMpWlWtDm2(:)  ! Avg wet wall fraction IC
        
        ! Work arrays from SAVE (task persistence)
        integer, allocatable :: ICpBtDm(:)     ! Bottom compartment per domain
        integer, allocatable :: ICpTpWaSrDm(:) ! Top compartment water storage per domain
        real(8), allocatable :: ArMpTpDm(:)    ! Area at top per domain
        real(8), allocatable :: AwlCorFac(:)   ! Correction factor
        real(8), allocatable :: FrMpWalWet(:,:)    ! Wet macropore wall fraction
        real(8), allocatable :: KDCrRlRef(:)   ! Reference crack conductivity
        real(8), allocatable :: QExcMtxDmCp(:,:)   ! Exchange flux per domain/comp
        real(8), allocatable :: QInIntSatDmCp(:,:) ! Interflow in per domain/comp
        real(8), allocatable :: QInMtxSatDmCp(:,:) ! Matrix sat inflow per domain/comp
        real(8), allocatable :: QInTopLatDm(:)     ! Top lateral inflow per domain
        real(8), allocatable :: QInTopVrtDm(:)     ! Top vertical inflow per domain
        real(8), allocatable :: QOutDrRapCp(:)     ! Rapid drainage out per comp
        real(8), allocatable :: QOutMtxSatDmCp(:,:)    ! Matrix sat outflow per domain/comp
        real(8), allocatable :: QOutMtxUnsDmCp(:,:)    ! Matrix unsat outflow per domain/comp
        real(8), allocatable :: SorpDmCp(:,:)      ! Sorptivity per domain/comp
        real(8), allocatable :: ThtSrpRefDmCp(:,:) ! Reference theta for sorption
        real(8), allocatable :: TimAbsCumDmCp(:,:) ! Cumulative absorption time
        real(8), allocatable :: VlMpDm(:)          ! Macropore volume per domain
        real(8), allocatable :: VlMpDmCp(:,:)      ! Macropore volume per domain/comp
        real(8), allocatable :: WaSrMpDm(:)        ! Water storage per domain
        real(8), allocatable :: WaSrMpDmCp(:,:)    ! Water storage per domain/comp
        real(8), allocatable :: ZBtDm(:)           ! Bottom depth per domain
        real(8), allocatable :: ZWaLevDm(:)        ! Water level per domain
        logical, allocatable :: flDraTub(:)        ! Drain tube flag per level
        logical, allocatable :: FlEndSrpEvt(:,:)   ! End sorption event flag
        
        ! State tracking
        real(8) :: WaSrMp = 0.0d0              ! Total water storage in macropores
        integer :: ICpBtPerZon = 0             ! Bottom compartment percolation zone
        integer :: ICpSatGWl = 0               ! Saturated compartment at GWL
        integer :: ICpSatPeGWl = 0             ! Saturated perched GWL compartment
        integer :: ICpTpPerZon = 0             ! Top compartment percolation zone
        integer :: ICpTpSatZon = 0             ! Top compartment saturated zone
        integer :: NnCrAr = 0                  ! Number of crack areas
        logical :: flBegin = .true.            ! Beginning of simulation flag
        
        ! Domain configuration
        integer :: NumDm = 0                   ! Number of domains
        integer :: NumSbDm = 0                 ! Subdomains in IC
        integer :: IcTopMP = 0                 ! Top compartment with macropores
        integer :: NumLevRapDra = 0            ! Number of rapid drainage levels
        real(8) :: Z_Tp = 0.0d0                ! Top depth of macropores
        real(8) :: Z_St = 0.0d0                ! Bottom of static macropores
        real(8) :: Z_Ic = 0.0d0                ! Bottom of IC domain
        real(8) :: Z_Ah = 0.0d0                ! Bottom of A-horizon
        real(8) :: ArMpTp = 0.0d0              ! Area fraction at top of macropores
        real(8) :: ArMpSs = 0.0d0              ! Area fraction at soil surface
        real(8) :: KsatCovLay = 0.0d0          ! Saturated K of covering layer
        real(8) :: KsMpSs = 0.0d0              ! Vertical K of macropores at surface
        real(8) :: PpIcTpMp = 0.0d0            ! Proportion IC at top macropores
        real(8) :: dtold = 0.0d0               ! Previous timestep length
        
        ! Groundwater tracking
        real(8) :: GWlFlCpZo = 0.0d0           ! GWL of full capillary zone
        integer :: NodGWlFlCpZo = 0
        real(8) :: ZDraBas = 0.0d0             ! Drainage basis level
        
        ! Iteration control
        integer :: IDecMpRat = 0               ! Convergence iteration counter
        logical :: FlDecMpRat = .false.        ! Decrease macropore fluxes
        logical :: flInitDraBas = .false.      ! Initialize drainage basis
        
        ! Flags
        logical :: flmacropore = .false.
    end type macropore_state_t

    ! ===========================================================================
    ! Solute Transport State
    ! ===========================================================================
    type :: solute_state_t
        ! Configuration switches
        integer :: swsolu = 0                  ! Switch for solute transport simulation
        integer :: swsp = 0                    ! Switch for sorption simulation
        integer :: swbr = 0                    ! Switch for breakthrough curve
        integer :: swbotbc = 1                 ! Switch for bottom BC
        integer :: nconc = 0                   ! Number of initial concentrations
        
        ! Concentrations
        real(8), allocatable :: cml(:)         ! Mobile concentration (mg/cm3)
        real(8), allocatable :: cmsy(:)        ! Total concentration
        real(8) :: cpond = 0.0d0               ! Ponding concentration
        real(8) :: csurf = 0.0d0               ! Surface solute amount
        real(8) :: cdrain = 0.0d0              ! Drainage concentration
        real(8) :: cseep = 0.0d0               ! Seepage concentration
        real(8) :: cpre = 0.0d0                ! Precipitation concentration
        real(8) :: cirr = 0.0d0                ! Irrigation concentration
        real(8) :: cref = 1.0d0                ! Reference concentration for Freundlich
        
        ! Cumulative amounts
        real(8) :: sampro = 0.0d0              ! Total in profile
        real(8) :: samini = 0.0d0              ! Initial in profile
        real(8) :: sqbot = 0.0d0               ! Through bottom
        real(8) :: sqdra = 0.0d0               ! To drainage
        real(8) :: sqprec = 0.0d0              ! In precipitation
        real(8) :: sqirrig = 0.0d0             ! In irrigation
        real(8) :: sqsur = 0.0d0               ! To surface water
        real(8) :: sqrap = 0.0d0               ! In rapid drainage
        real(8) :: dectot = 0.0d0              ! Decomposition
        real(8) :: rottot = 0.0d0              ! Root extraction
        real(8) :: solbal = 0.0d0              ! Balance error
        
        ! Intermediate amounts
        real(8) :: imsqbot = 0.0d0
        real(8) :: imsqdra = 0.0d0
        real(8) :: imsqprec = 0.0d0
        real(8) :: imsqirrig = 0.0d0
        real(8) :: imdectot = 0.0d0
        real(8) :: imrottot = 0.0d0
        real(8) :: isqbot = 0.0d0
        real(8) :: isqtop = 0.0d0
        
        ! Macropore solute
        real(8) :: samcra = 0.0d0              ! In cracks
        
        ! Age tracer
        real(8) :: AgeGwl1m = 0.0d0            ! GW age upper 1m sat zone
        real(8) :: icAgeBot = 0.0d0
        real(8) :: icAgeRot = 0.0d0
        real(8) :: icAgeSur = 0.0d0
        real(8), allocatable :: icAgeDra(:)
        
        ! Age tracer boundary/pond state (persists between timesteps)
        real(8) :: Ageirr = 0.0d0              ! Age of irrigation water
        real(8) :: Agedrain = 0.0d0            ! Age of drainage water
        real(8) :: Agepre = 0.0d0              ! Age of precipitation
        real(8) :: Agepond = 0.0d0             ! Age of ponding water
        real(8) :: Agepondm1 = 0.0d0           ! Age of ponding (prev timestep)
        real(8) :: icAgetopupw = 0.0d0         ! Incremental age leaving top
        real(8) :: icAgetopdwn = 0.0d0         ! Incremental age entering top
        real(8) :: ArMpSs = 0.0d0              ! Area fraction macropores at surface
        
        ! Transport parameters
        real(8) :: ddif = 0.0d0                ! Molecular diffusion coefficient
        real(8) :: frexp = 1.0d0               ! Freundlich exponent
        real(8) :: tscf = 1.0d0                ! Relative uptake by roots
        real(8) :: dtsolu = 0.0d0              ! Max time step for solute
        
        ! Decomposition parameters
        real(8) :: gampar = 0.0d0              ! Temp reduction factor
        real(8) :: bexp = 0.7d0                ! Dryness exponent
        real(8) :: rtheta = 0.01d0             ! Min theta for decomposition
        real(8) :: decsat = 0.0d0              ! Decomposition in aquifer
        
        ! Aquifer parameters for breakthrough
        real(8) :: daquif = 0.0d0              ! Aquifer thickness
        real(8) :: poros = 0.0d0               ! Aquifer porosity
        real(8) :: kfsat = 0.0d0               ! Adsorption in aquifer
        
        ! Salt stress parameters
        real(8) :: salthead = 0.0d0            ! Salt to osmotic head
        real(8) :: saltmax = 0.0d0             ! Threshold concentration
        real(8) :: saltslope = 0.0d0           ! Uptake decline
        
        ! Parameters per layer
        real(8), allocatable :: ldis(:)        ! Dispersion length
        real(8), allocatable :: kf(:)          ! Freundlich coefficient
        real(8), allocatable :: decpot(:)      ! Potential decomposition rate
        real(8), allocatable :: fdepth(:)      ! Depth reduction factor
        
        ! Tables
        real(8), allocatable :: cseeptab(:)    ! Seepage concentration table
        real(8), allocatable :: zc(:)          ! Depths for initial concentrations
        
        ! Flags
        logical :: flsolute = .false.
        logical :: flAgeTracer = .false.
    end type solute_state_t

    ! ===========================================================================
    ! Heat Flow State
    ! ===========================================================================
    type :: heat_state_t
        ! Configuration switches
        integer :: swhea = 0                   ! Switch for heat flow simulation
        integer :: swcalt = 1                  ! Method: 1=analytical, 2=numerical
        integer :: swtopbhea = 1               ! Top BC: 1=air temp, 2=measured
        integer :: swbotbhea = 1               ! Bottom BC: 1=zero flux, 2=prescribed
        integer :: swfrost = 0                 ! Switch for frost reduction
        integer :: nheat = 0                   ! Number of initial temperatures
        
        ! Soil temperatures
        real(8), allocatable :: tsoil(:)       ! Temperature per compartment (°C)
        real(8) :: tetop = 0.0d0               ! Top temperature
        real(8) :: tebot = 0.0d0               ! Bottom temperature
        
        ! Thermal properties per compartment
        real(8), allocatable :: heacap(:)      ! Heat capacity (J/cm3/K)
        real(8), allocatable :: heacon(:)      ! Heat conductivity (J/cm/K/d)
        
        ! Frost reduction factor per compartment
        real(8), allocatable :: rfcp(:)        ! Reduction factor for frozen conditions
        
        ! Soil composition per compartment
        real(8), allocatable :: fclay(:)       ! Clay content
        real(8), allocatable :: forg(:)        ! Organic matter content
        real(8), allocatable :: fquartz(:)     ! Sand+silt content
        
        ! Soil composition per layer
        real(8), allocatable :: pclay(:)       ! Clay content per layer
        real(8), allocatable :: psand(:)       ! Sand content per layer
        real(8), allocatable :: psilt(:)       ! Silt content per layer
        real(8), allocatable :: orgmat(:)      ! Organic matter per layer
        
        ! Boundary conditions
        real(8) :: tmean = 0.0d0               ! Mean annual surface temperature
        real(8) :: tampli = 0.0d0              ! Surface temperature amplitude
        real(8) :: timref = 0.0d0              ! Time of max temperature
        real(8) :: ddamp = 0.0d0               ! Damping depth
        
        ! Boundary condition tables
        real(8), allocatable :: tembtab(:)     ! Bottom temperature table
        real(8), allocatable :: temtoptab(:)   ! Top temperature table
        real(8), allocatable :: zh(:)          ! Depths for initial temperatures
        
        ! Frost
        real(8) :: zfrosttop = 0.0d0           ! Top of frost layer
        real(8) :: zfrostbot = 0.0d0           ! Bottom of frost layer
        real(8) :: tfroststa = 0.0d0           ! Frost start temperature
        real(8) :: tfrostend = 0.0d0           ! Frost end temperature
        integer :: nodfrostbot = 0             ! Deepest frost node
        
        ! Flags
        logical :: fltemperature = .false.
    end type heat_state_t

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
    ! Surface Water State (for extended drainage)
    ! ===========================================================================
    type :: surfacewater_state_t
        ! Water levels
        real(8) :: wlp = 0.0d0                 ! Primary water level
        real(8) :: wls = 0.0d0                 ! Secondary water level
        real(8) :: wlsold = 0.0d0              ! Previous sec water level
        real(8) :: wlstar = 0.0d0              ! Target water level
        real(8) :: hwlman = 0.0d0              ! Managed water level
        real(8) :: vtair = 0.0d0               ! Total air volume in profile
        real(8), allocatable :: wlsbak(:)      ! Last 4 water levels for oscillation check
        
        ! Storage
        real(8) :: swst = 0.0d0                ! Surface water storage
        real(8) :: swstini = 0.0d0             ! Initial storage
        
        ! Fluxes
        real(8) :: qdrd = 0.0d0                ! Drainage discharge to secondary
        real(8) :: cqdrd = 0.0d0               ! Cumulative drainage discharge
        real(8) :: cwsupp = 0.0d0              ! Cumulative water supply
        real(8) :: cwout = 0.0d0               ! Cumulative water out
        real(8) :: runots = 0.0d0              ! Runoff during timestep
        real(8) :: QRapDra = 0.0d0             ! Rapid drainage flux
        
        ! Management parameters
        integer :: imper = 0                   ! Current management period
        integer :: nmper = 0                   ! Number of management periods
        integer :: numadj = 0                  ! Number of adjustments
        integer :: swsrf = 0                   ! Surface water switch
        integer :: swsec = 0                   ! Secondary water switch
        integer :: swqhr = 0                   ! Q-h relation switch
        real(8) :: osswlm = 0.0d0              ! Oscillation tolerance
        real(8), allocatable :: impend(:)      ! Management period end dates
        integer, allocatable :: swman(:)       ! Management type per period
        real(8), allocatable :: hbweir(:)      ! Weir crest level per period
        real(8), allocatable :: wldip(:)       ! Dip below target for supply
        real(8), allocatable :: alphaw(:)      ! Weir discharge coef a
        real(8), allocatable :: betaw(:)       ! Weir discharge coef b
        real(8), allocatable :: wscap(:)       ! Max supply capacity per period
        real(8), allocatable :: dropr(:)       ! Max drop rate per period
        real(8), allocatable :: intwl(:)       ! Adjustment interval per period
        integer, allocatable :: nphase(:)      ! Number of phases per period
        real(8), allocatable :: gwlcrit(:,:)   ! Critical gwl per period/phase
        real(8), allocatable :: wlsman(:,:)    ! Target level per period/phase
        real(8), allocatable :: vcrit(:,:)     ! Critical air volume per period/phase
        real(8), allocatable :: hcrit(:,:)     ! Critical pressure head per period/phase
        integer, allocatable :: nodhd(:)       ! Node for head criterion per period
        
        ! Lookup tables
        real(8), allocatable :: wlstab(:)      ! Water level vs time table
        real(8), allocatable :: wlptab(:)      ! Primary water level table
        real(8), allocatable :: sttab(:,:)     ! Storage-level table
        real(8), allocatable :: owltab(:,:)    ! Open water level tables per drain level
        real(8), allocatable :: qqhtab(:,:)    ! Q-H table per period
        
        ! Flags
        logical :: flsurfacewater = .false.
        logical :: overfl = .false.
        logical :: fldecdt = .false.           ! Decrease timestep flag
        logical :: fldtmin = .false.           ! At minimum timestep flag
    end type surfacewater_state_t

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
    
    subroutine soil_state_init(soil, numnod, numlay)
        type(soil_state_t), intent(inout) :: soil
        integer, intent(in) :: numnod, numlay
        
        integer :: alloc_stat
        
        ! Allocate node-based arrays
        allocate(soil%h(numnod), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call log_error('soil_init', 'Failed to allocate h array')
            return
        end if
        allocate(soil%hm1(numnod))
        allocate(soil%theta(numnod))
        allocate(soil%thetm1(numnod))
        allocate(soil%k(numnod+1))
        allocate(soil%kmean(numnod+1))
        allocate(soil%dimoca(numnod))
        allocate(soil%z(numnod))
        allocate(soil%dz(numnod))
        allocate(soil%ztopcp(numnod))
        allocate(soil%zbotcp(numnod))
        allocate(soil%disnod(numnod+1))
        allocate(soil%thetar(numnod))
        allocate(soil%thetas(numnod))
        allocate(soil%hroot(numnod))
        allocate(soil%twilt(numnod))
        allocate(soil%rfcp(numnod))
        allocate(soil%q(numnod+1))
        allocate(soil%evp(numnod))
        allocate(soil%qrot(numnod))
        allocate(soil%qpotrot(numnod))
        allocate(soil%qssdi(numnod))
        allocate(soil%qredwet(numnod))
        allocate(soil%qreddry(numnod))
        allocate(soil%qredsol(numnod))
        allocate(soil%qredfrs(numnod))
        allocate(soil%inq(numnod+1))
        allocate(soil%inqrot(numnod))
        allocate(soil%inqssdi(numnod))
        allocate(soil%ithetabeg(numnod))
        allocate(soil%layer(numnod))
        allocate(soil%indeks(numnod))
        
        ! Allocate layer-based arrays
        allocate(soil%bdens(numlay))
        allocate(soil%cofani(numlay))
        allocate(soil%ksatfit(numlay))
        allocate(soil%ksatexm(numlay))
        allocate(soil%thetsl(numlay))
        allocate(soil%botcom(numlay))
        allocate(soil%nod1lay(numlay))
        allocate(soil%paramvg(21, numlay))
        allocate(soil%cofgen(21, numnod))
        
        
        ! Initialize to zero
        soil%h = 0.0d0
        soil%hm1 = 0.0d0
        soil%theta = 0.0d0
        soil%thetm1 = 0.0d0
        soil%k = 0.0d0
        soil%kmean = 0.0d0
        soil%dimoca = 0.0d0
        soil%z = 0.0d0
        soil%dz = 0.0d0
        soil%numnod = numnod
        soil%numlay = numlay
        
    end subroutine soil_state_init
    
    subroutine drainage_state_init(drain, nrlevs, numnod)
        type(drainage_state_t), intent(inout) :: drain
        integer, intent(in) :: nrlevs, numnod
        
        ! Drainage fluxes per level
        allocate(drain%qdrain(nrlevs))
        allocate(drain%cqdrain(nrlevs))
        allocate(drain%cqdrainin(nrlevs))
        allocate(drain%cqdrainout(nrlevs))
        allocate(drain%drainl(nrlevs))
        
        ! Spatial arrays
        allocate(drain%qdra(nrlevs, numnod))
        allocate(drain%inqdra(nrlevs, numnod))
        allocate(drain%inqdra_in(nrlevs, numnod))
        allocate(drain%inqdra_out(nrlevs, numnod))
        allocate(drain%qdraincomp(numnod))
        
        ! Resistances and geometry per level
        allocate(drain%drares(nrlevs))
        allocate(drain%infres(nrlevs))
        allocate(drain%L(nrlevs))
        allocate(drain%wetper(nrlevs))
        allocate(drain%zbotdr(nrlevs))
        allocate(drain%rdrain(nrlevs))
        allocate(drain%rinfi(nrlevs))
        allocate(drain%rentry(nrlevs))
        allocate(drain%rexit(nrlevs))
        allocate(drain%gwlinf(nrlevs))
        allocate(drain%widthr(nrlevs))
        allocate(drain%taludr(nrlevs))
        allocate(drain%swallo(nrlevs))
        allocate(drain%swdtyp(nrlevs))
        allocate(drain%swtopdislay(nrlevs))
        allocate(drain%zTopDisLay(nrlevs))
        allocate(drain%fTopDisLay(nrlevs))
        
        ! Initialize to zero
        drain%qdrain = 0.0d0
        drain%cqdrain = 0.0d0
        drain%cqdrainin = 0.0d0
        drain%cqdrainout = 0.0d0
        drain%drainl = 0.0d0
        drain%qdra = 0.0d0
        drain%inqdra = 0.0d0
        drain%inqdra_in = 0.0d0
        drain%inqdra_out = 0.0d0
        drain%qdraincomp = 0.0d0
        drain%drares = 0.0d0
        drain%infres = 0.0d0
        drain%L = 0.0d0
        drain%wetper = 0.0d0
        drain%zbotdr = 0.0d0
        drain%rdrain = 0.0d0
        drain%rinfi = 0.0d0
        drain%rentry = 0.0d0
        drain%rexit = 0.0d0
        drain%gwlinf = 0.0d0
        drain%widthr = 0.0d0
        drain%taludr = 0.0d0
        drain%swallo = 0
        drain%swdtyp = 0
        drain%swtopdislay = 0
        drain%zTopDisLay = 0.0d0
        drain%fTopDisLay = 0.0d0
        drain%nrlevs = nrlevs
        
    end subroutine drainage_state_init
    
    subroutine surfacewater_state_init(swstate, nmper, mamte, maowl, nrlevs)
        type(surfacewater_state_t), intent(inout) :: swstate
        integer, intent(in) :: nmper   ! Max management periods
        integer, intent(in) :: mamte   ! Max meteo entries for gwl phases
        integer, intent(in) :: maowl   ! Max open water level entries
        integer, intent(in) :: nrlevs  ! Number of drainage levels
        
        ! Water level history for oscillation check
        allocate(swstate%wlsbak(4))
        
        ! Management period arrays
        allocate(swstate%impend(nmper))
        allocate(swstate%swman(nmper))
        allocate(swstate%hbweir(nmper))
        allocate(swstate%wldip(nmper))
        allocate(swstate%alphaw(nmper))
        allocate(swstate%betaw(nmper))
        allocate(swstate%wscap(nmper))
        allocate(swstate%dropr(nmper))
        allocate(swstate%intwl(nmper))
        allocate(swstate%nphase(nmper))
        allocate(swstate%nodhd(nmper))
        allocate(swstate%gwlcrit(nmper, mamte))
        allocate(swstate%wlsman(nmper, mamte))
        allocate(swstate%vcrit(nmper, mamte))
        allocate(swstate%hcrit(nmper, mamte))
        
        ! Lookup tables
        allocate(swstate%wlstab(2*maowl))
        allocate(swstate%wlptab(2*maowl))
        allocate(swstate%sttab(22, 2))
        allocate(swstate%owltab(nrlevs, 2*maowl))
        allocate(swstate%qqhtab(nmper, 22))
        
        ! Initialize
        swstate%wlsbak = 0.0d0
        swstate%impend = 0.0d0
        swstate%swman = 0
        swstate%hbweir = 0.0d0
        swstate%wldip = 0.0d0
        swstate%alphaw = 0.0d0
        swstate%betaw = 0.0d0
        swstate%wscap = 0.0d0
        swstate%dropr = 0.0d0
        swstate%intwl = 0.0d0
        swstate%nphase = 0
        swstate%nodhd = 0
        swstate%gwlcrit = 0.0d0
        swstate%wlsman = 0.0d0
        swstate%vcrit = 0.0d0
        swstate%hcrit = 0.0d0
        swstate%wlstab = 0.0d0
        swstate%wlptab = 0.0d0
        swstate%sttab = 0.0d0
        swstate%owltab = 0.0d0
        swstate%qqhtab = 0.0d0
        swstate%nmper = nmper
        
    end subroutine surfacewater_state_init
    
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
    
    subroutine macropore_state_init(macro, numnod, maxdom, maxdra)
        type(macropore_state_t), intent(inout) :: macro
        integer, intent(in) :: numnod, maxdom, maxdra
        
        ! Original per-compartment arrays
        allocate(macro%DiPoCp(numnod))
        allocate(macro%FrArMtrx(numnod))
        allocate(macro%VlMpDyCp(numnod))
        allocate(macro%VlMpStCp(numnod))
        allocate(macro%VlMpStDm1(numnod))
        allocate(macro%VlMpStDm2(numnod))
        allocate(macro%QExcMpMtx(numnod))
        allocate(macro%SubsidCp(numnod))
        allocate(macro%PpDmCp(maxdom, numnod))
        allocate(macro%dFdhMp(numnod))
        allocate(macro%iQOutDrRapCp(numnod))
        allocate(macro%iQExcMtxDm1Cp(numnod))
        allocate(macro%iQExcMtxDm2Cp(numnod))
        allocate(macro%IAvFrMpWlWtDm1(numnod))
        allocate(macro%IAvFrMpWlWtDm2(numnod))
        
        ! Work arrays from SAVE (task persistence)
        allocate(macro%ICpBtDm(maxdom))
        allocate(macro%ICpTpWaSrDm(maxdom))
        allocate(macro%ArMpTpDm(maxdom))
        allocate(macro%AwlCorFac(numnod))
        allocate(macro%FrMpWalWet(maxdom, numnod))
        allocate(macro%KDCrRlRef(maxdra))
        allocate(macro%QExcMtxDmCp(maxdom, numnod))
        allocate(macro%QInIntSatDmCp(maxdom, numnod))
        allocate(macro%QInMtxSatDmCp(maxdom, numnod))
        allocate(macro%QInTopLatDm(maxdom))
        allocate(macro%QInTopVrtDm(maxdom))
        allocate(macro%QOutDrRapCp(numnod))
        allocate(macro%QOutMtxSatDmCp(maxdom, numnod))
        allocate(macro%QOutMtxUnsDmCp(maxdom, numnod))
        allocate(macro%SorpDmCp(maxdom, numnod))
        allocate(macro%ThtSrpRefDmCp(maxdom, numnod))
        allocate(macro%TimAbsCumDmCp(maxdom, numnod))
        allocate(macro%VlMpDm(maxdom))
        allocate(macro%VlMpDmCp(maxdom, numnod))
        allocate(macro%WaSrMpDm(maxdom))
        allocate(macro%WaSrMpDmCp(maxdom, numnod))
        allocate(macro%ZBtDm(maxdom))
        allocate(macro%ZWaLevDm(maxdom))
        allocate(macro%flDraTub(maxdra))
        allocate(macro%FlEndSrpEvt(maxdom, numnod))
        
        ! Initialize per-compartment arrays
        macro%DiPoCp = 0.0d0
        macro%FrArMtrx = 1.0d0
        macro%VlMpDyCp = 0.0d0
        macro%VlMpStCp = 0.0d0
        macro%VlMpStDm1 = 0.0d0
        macro%VlMpStDm2 = 0.0d0
        macro%QExcMpMtx = 0.0d0
        macro%dFdhMp = 0.0d0
        macro%iQOutDrRapCp = 0.0d0
        macro%iQExcMtxDm1Cp = 0.0d0
        macro%iQExcMtxDm2Cp = 0.0d0
        macro%IAvFrMpWlWtDm1 = 0.0d0
        macro%IAvFrMpWlWtDm2 = 0.0d0
        
        ! Initialize work arrays
        macro%ICpBtDm = 0
        macro%ICpTpWaSrDm = 0
        macro%ArMpTpDm = 0.0d0
        macro%AwlCorFac = 0.0d0
        macro%FrMpWalWet = 0.0d0
        macro%KDCrRlRef = 0.0d0
        macro%QExcMtxDmCp = 0.0d0
        macro%QInIntSatDmCp = 0.0d0
        macro%QInMtxSatDmCp = 0.0d0
        macro%QInTopLatDm = 0.0d0
        macro%QInTopVrtDm = 0.0d0
        macro%QOutDrRapCp = 0.0d0
        macro%QOutMtxSatDmCp = 0.0d0
        macro%QOutMtxUnsDmCp = 0.0d0
        macro%SorpDmCp = 0.0d0
        macro%ThtSrpRefDmCp = 0.0d0
        macro%TimAbsCumDmCp = 0.0d0
        macro%VlMpDm = 0.0d0
        macro%VlMpDmCp = 0.0d0
        macro%WaSrMpDm = 0.0d0
        macro%WaSrMpDmCp = 0.0d0
        macro%ZBtDm = 0.0d0
        macro%ZWaLevDm = 0.0d0
        macro%flDraTub = .false.
        macro%FlEndSrpEvt = .false.
        
    end subroutine macropore_state_init
    
    subroutine solute_state_init(solute, numnod, numlay, nrlevs)
        type(solute_state_t), intent(inout) :: solute
        integer, intent(in) :: numnod, numlay, nrlevs
        
        allocate(solute%cml(numnod))
        allocate(solute%cmsy(numnod))
        allocate(solute%ldis(numlay))
        allocate(solute%kf(numlay))
        allocate(solute%decpot(numlay))
        allocate(solute%fdepth(numlay))
        allocate(solute%icAgeDra(nrlevs))
        allocate(solute%cseeptab(2*MABBC))
        allocate(solute%zc(numnod))
        
        solute%cml = 0.0d0
        solute%cmsy = 0.0d0
        solute%ldis = 0.0d0
        solute%kf = 0.0d0
        solute%decpot = 0.0d0
        solute%fdepth = 1.0d0
        solute%icAgeDra = 0.0d0
        solute%cseeptab = 0.0d0
        solute%zc = 0.0d0
        
        
    end subroutine solute_state_init
    
    subroutine heat_state_init(heat, numnod, numlay)
        type(heat_state_t), intent(inout) :: heat
        integer, intent(in) :: numnod, numlay
        
        allocate(heat%tsoil(numnod))
        allocate(heat%heacap(numnod))
        allocate(heat%heacon(numnod))
        allocate(heat%rfcp(numnod))
        allocate(heat%fclay(numnod))
        allocate(heat%forg(numnod))
        allocate(heat%fquartz(numnod))
        allocate(heat%pclay(numlay))
        allocate(heat%psand(numlay))
        allocate(heat%psilt(numlay))
        allocate(heat%orgmat(numlay))
        allocate(heat%tembtab(2*MABBC))
        allocate(heat%temtoptab(2*MABBC))
        allocate(heat%zh(numnod))
        
        heat%tsoil = 10.0d0  ! Default temperature
        heat%heacap = 0.0d0
        heat%heacon = 0.0d0
        heat%rfcp = 1.0d0    ! No reduction by default
        heat%fclay = 0.0d0
        heat%forg = 0.0d0
        heat%fquartz = 0.0d0
        heat%tembtab = 0.0d0
        heat%temtoptab = 0.0d0
        heat%zh = 0.0d0
        
        
    end subroutine heat_state_init
    
    subroutine boundary_state_init(boundary, nday)
        type(boundary_state_t), intent(inout) :: boundary
        integer, intent(in), optional :: nday
        
        integer :: nd
        
        ! Default number of days
        nd = MADAY
        if (present(nday)) nd = nday
        
        ! Allocate tables with reasonable default sizes
        ! These may need to be reallocated based on input
        allocate(boundary%gwltab(2*MABBC))
        allocate(boundary%haqtab(2*MABBC))
        allocate(boundary%qbotab(2*MABBC))
        allocate(boundary%hbotab(2*MABBC))
        allocate(boundary%runonarr(nd))
        allocate(boundary%pondmxtab(2*MAIRG))
        
        boundary%gwltab = 0.0d0
        boundary%haqtab = 0.0d0
        boundary%qbotab = 0.0d0
        boundary%hbotab = 0.0d0
        boundary%runonarr = 0.0d0
        boundary%pondmxtab = 0.0d0
        
    end subroutine boundary_state_init
    
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
    
    !> @brief Finalize soil state (deallocate arrays)
    subroutine soil_state_finalize(soil)
        type(soil_state_t), intent(inout) :: soil
        
        ! Node-based arrays
        if (allocated(soil%h)) deallocate(soil%h)
        if (allocated(soil%theta)) deallocate(soil%theta)
        if (allocated(soil%k)) deallocate(soil%k)
        if (allocated(soil%z)) deallocate(soil%z)
        if (allocated(soil%dz)) deallocate(soil%dz)
        if (allocated(soil%hm1)) deallocate(soil%hm1)
        if (allocated(soil%thetm1)) deallocate(soil%thetm1)
        if (allocated(soil%q)) deallocate(soil%q)
        if (allocated(soil%qrot)) deallocate(soil%qrot)
        if (allocated(soil%kmean)) deallocate(soil%kmean)
        if (allocated(soil%rfcp)) deallocate(soil%rfcp)
        
        ! Layer-based arrays
        if (allocated(soil%bdens)) deallocate(soil%bdens)
        if (allocated(soil%paramvg)) deallocate(soil%paramvg)
        
    end subroutine soil_state_finalize
    
    !> @brief Finalize drainage state (deallocate arrays)
    subroutine drain_state_finalize(drain)
        type(drainage_state_t), intent(inout) :: drain
        
        ! Flux arrays
        if (allocated(drain%qdrain)) deallocate(drain%qdrain)
        if (allocated(drain%cqdrain)) deallocate(drain%cqdrain)
        if (allocated(drain%cqdrainin)) deallocate(drain%cqdrainin)
        if (allocated(drain%cqdrainout)) deallocate(drain%cqdrainout)
        if (allocated(drain%drainl)) deallocate(drain%drainl)
        
        ! Spatial arrays
        if (allocated(drain%qdra)) deallocate(drain%qdra)
        if (allocated(drain%inqdra)) deallocate(drain%inqdra)
        if (allocated(drain%inqdra_in)) deallocate(drain%inqdra_in)
        if (allocated(drain%inqdra_out)) deallocate(drain%inqdra_out)
        if (allocated(drain%qdraincomp)) deallocate(drain%qdraincomp)
        
        ! Resistance/geometry arrays
        if (allocated(drain%drares)) deallocate(drain%drares)
        if (allocated(drain%infres)) deallocate(drain%infres)
        if (allocated(drain%L)) deallocate(drain%L)
        if (allocated(drain%wetper)) deallocate(drain%wetper)
        if (allocated(drain%zbotdr)) deallocate(drain%zbotdr)
        if (allocated(drain%rdrain)) deallocate(drain%rdrain)
        if (allocated(drain%rinfi)) deallocate(drain%rinfi)
        if (allocated(drain%rentry)) deallocate(drain%rentry)
        if (allocated(drain%rexit)) deallocate(drain%rexit)
        if (allocated(drain%gwlinf)) deallocate(drain%gwlinf)
        if (allocated(drain%widthr)) deallocate(drain%widthr)
        if (allocated(drain%taludr)) deallocate(drain%taludr)
        if (allocated(drain%swallo)) deallocate(drain%swallo)
        if (allocated(drain%swdtyp)) deallocate(drain%swdtyp)
        if (allocated(drain%swtopdislay)) deallocate(drain%swtopdislay)
        if (allocated(drain%zTopDisLay)) deallocate(drain%zTopDisLay)
        if (allocated(drain%fTopDisLay)) deallocate(drain%fTopDisLay)
        
    end subroutine drain_state_finalize
    
    !> @brief Finalize surfacewater state (deallocate arrays)
    subroutine surfacewater_state_finalize(swstate)
        type(surfacewater_state_t), intent(inout) :: swstate
        
        if (allocated(swstate%wlsbak)) deallocate(swstate%wlsbak)
        if (allocated(swstate%impend)) deallocate(swstate%impend)
        if (allocated(swstate%swman)) deallocate(swstate%swman)
        if (allocated(swstate%hbweir)) deallocate(swstate%hbweir)
        if (allocated(swstate%wldip)) deallocate(swstate%wldip)
        if (allocated(swstate%alphaw)) deallocate(swstate%alphaw)
        if (allocated(swstate%betaw)) deallocate(swstate%betaw)
        if (allocated(swstate%wscap)) deallocate(swstate%wscap)
        if (allocated(swstate%dropr)) deallocate(swstate%dropr)
        if (allocated(swstate%intwl)) deallocate(swstate%intwl)
        if (allocated(swstate%nphase)) deallocate(swstate%nphase)
        if (allocated(swstate%nodhd)) deallocate(swstate%nodhd)
        if (allocated(swstate%gwlcrit)) deallocate(swstate%gwlcrit)
        if (allocated(swstate%wlsman)) deallocate(swstate%wlsman)
        if (allocated(swstate%vcrit)) deallocate(swstate%vcrit)
        if (allocated(swstate%hcrit)) deallocate(swstate%hcrit)
        if (allocated(swstate%wlstab)) deallocate(swstate%wlstab)
        if (allocated(swstate%wlptab)) deallocate(swstate%wlptab)
        if (allocated(swstate%sttab)) deallocate(swstate%sttab)
        if (allocated(swstate%owltab)) deallocate(swstate%owltab)
        if (allocated(swstate%qqhtab)) deallocate(swstate%qqhtab)
        
    end subroutine surfacewater_state_finalize
    
    !> @brief Finalize boundary state (deallocate arrays)
    subroutine boundary_state_finalize(boundary)
        type(boundary_state_t), intent(inout) :: boundary
        
        if (allocated(boundary%gwltab)) deallocate(boundary%gwltab)
        if (allocated(boundary%haqtab)) deallocate(boundary%haqtab)
        if (allocated(boundary%qbotab)) deallocate(boundary%qbotab)
        if (allocated(boundary%hbotab)) deallocate(boundary%hbotab)
        if (allocated(boundary%runonarr)) deallocate(boundary%runonarr)
        if (allocated(boundary%pondmxtab)) deallocate(boundary%pondmxtab)
        
    end subroutine boundary_state_finalize
    
    !> @brief Finalize solute state (deallocate arrays)
    subroutine solute_state_finalize(solute)
        type(solute_state_t), intent(inout) :: solute
        
        if (allocated(solute%cml)) deallocate(solute%cml)
        if (allocated(solute%cmsy)) deallocate(solute%cmsy)
        if (allocated(solute%icAgeDra)) deallocate(solute%icAgeDra)
        if (allocated(solute%ldis)) deallocate(solute%ldis)
        if (allocated(solute%kf)) deallocate(solute%kf)
        if (allocated(solute%decpot)) deallocate(solute%decpot)
        if (allocated(solute%fdepth)) deallocate(solute%fdepth)
        if (allocated(solute%cseeptab)) deallocate(solute%cseeptab)
        if (allocated(solute%zc)) deallocate(solute%zc)
        
    end subroutine solute_state_finalize
    
    !> @brief Finalize heat state (deallocate arrays)
    subroutine heat_state_finalize(heat)
        type(heat_state_t), intent(inout) :: heat
        
        if (allocated(heat%tsoil)) deallocate(heat%tsoil)
        if (allocated(heat%heacap)) deallocate(heat%heacap)
        if (allocated(heat%heacon)) deallocate(heat%heacon)
        if (allocated(heat%rfcp)) deallocate(heat%rfcp)
        if (allocated(heat%fclay)) deallocate(heat%fclay)
        if (allocated(heat%forg)) deallocate(heat%forg)
        if (allocated(heat%fquartz)) deallocate(heat%fquartz)
        if (allocated(heat%pclay)) deallocate(heat%pclay)
        if (allocated(heat%psand)) deallocate(heat%psand)
        if (allocated(heat%psilt)) deallocate(heat%psilt)
        if (allocated(heat%orgmat)) deallocate(heat%orgmat)
        if (allocated(heat%tembtab)) deallocate(heat%tembtab)
        if (allocated(heat%temtoptab)) deallocate(heat%temtoptab)
        if (allocated(heat%zh)) deallocate(heat%zh)
        
    end subroutine heat_state_finalize
    
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
        call drain_state_finalize(state%drain)
        call surfacewater_state_finalize(state%surfwater)
        call boundary_state_finalize(state%boundary)
        call solute_state_finalize(state%solute)
        call heat_state_finalize(state%heat)
        call irrigation_state_finalize(state%irrig)
        call tillage_state_finalize(state%tillage)
        call wofost_soil_state_finalize(state%wofost_soil)
        call oxygenstress_state_finalize(state%oxystress)
        
        state%initialized = .false.
        
    end subroutine swap_state_finalize

end module swap_state_mod
