!> Atmosphere State Module
!!
!! This module defines the atmosphere/meteorology state for SWAP.
!! Manages meteorological forcing data, evapotranspiration calculations,
!! precipitation, and interception processes.
!!
!! Key components:
!! - Configuration switches (swetr, swinter, etc.)
!! - Current meteorological conditions (temperature, radiation, etc.)
!! - Precipitation and interception state
!! - ET demand and actual fluxes
!! - Cumulative and intermediate flux tracking
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module atmosphere_state_mod
   implicit none
   private

    public :: atmosphere_state_t
    public :: atmosphere_state_init
    public :: atmosphere_state_finalize
    public :: atmosphere_state_reset_cumulative
    public :: atmosphere_state_reset_intermediate

    type :: atmosphere_state_t
        ! Meteo input configuration (loaded from .swp)
        integer :: swetr = 0                  !! 0=daily basic data, 1=ETref values
        integer :: swdivide = 0               !! Partitioning switch
        integer :: swmetdetail = 0            !! 0=daily, 1=detailed
        integer :: swmeteo = 0                !! 1=no detailed needed, 2=detailed needed
        integer :: swrain = 0                 !! Rain input mode
        integer :: swetsine = 0               !! ETSine switch
        integer :: swinter = 0                !! Interception method switch
        integer :: swMetFilAll = 0            !! Meteo-in-one-file switch
        integer :: nmetdetail = 0             !! Number of detailed records per day
        real(8) :: altw = 0.0d0               !! Wind measurement height (m)
        real(8) :: angstroma = 0.0d0          !! Angstrom coefficient a
        real(8) :: angstromb = 0.0d0          !! Angstrom coefficient b

        ! Current meteorological values
        real(8) :: tav = 0.0d0                 !! Average air temperature (°C)
        real(8) :: tavd = 0.0d0                !! Average daytime temperature
        real(8) :: tmn = 0.0d0                 !! Minimum temperature
        real(8) :: tmx = 0.0d0                 !! Maximum temperature
        real(8) :: tmnr = 0.0d0                !! 7-day average Tmin
        real(8) :: rad = 0.0d0                 !! Global solar radiation (J/m2/d)
        real(8) :: rh = 0.0d0                  !! Relative humidity (-)
        real(8) :: lat = 0.0d0                 !! Latitude (degrees)
        real(8) :: alt = 0.0d0                 !! Altitude (m)
        real(8) :: daylp = 0.0d0               !! Photoperiodic daylength (hours)
        
        ! Precipitation
        real(8) :: grai = 0.0d0                !! Gross rain flux (cm/d)
        real(8) :: graidt = 0.0d0              !! Gross precip during timestep
        real(8) :: nraida = 0.0d0              !! Net rain daily average
        real(8) :: nraidt = 0.0d0              !! Net rain during timestep
        real(8) :: finterception = 0.0d0       !! Net/gross rain ratio
        
        ! Evapotranspiration
        real(8) :: peva = 0.0d0                !! Potential soil evaporation (cm/d)
        real(8) :: pevaday = 0.0d0             !! Potential E of one day
        real(8) :: ptra = 0.0d0                !! Potential transpiration (cm/d)
        real(8) :: ptraday = 0.0d0             !! Potential T of one day
        real(8) :: tra = 0.0d0                 !! Actual transpiration (cm/d)
        real(8) :: reva = 0.0d0                !! Actual soil evaporation
        real(8) :: atmdem = 0.0d0              !! Atmospheric demand
        real(8) :: es0 = 0.0d0                 !! Potential E wet bare soil
        real(8) :: et0 = 0.0d0                 !! Potential T dry crop
        real(8) :: ew0 = 0.0d0                 !! Potential T wet crop
        
        ! Interception
        real(8) :: aintcdt = 0.0d0             !! Interception during timestep
        real(8) :: sicact = 0.0d0              !! Water stored on canopy
        real(8) :: siccapact = 0.0d0           !! Interception capacity
        
        ! Cumulative values
        real(8) :: cgrai = 0.0d0               !! Cumulative gross precipitation
        real(8) :: cnrai = 0.0d0               !! Cumulative net precipitation
        real(8) :: caintc = 0.0d0              !! Cumulative interception
        real(8) :: cevap = 0.0d0               !! Cumulative actual evaporation
        real(8) :: cpeva = 0.0d0               !! Cumulative potential evaporation
        real(8) :: cptra = 0.0d0               !! Cumulative potential transpiration
        
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

contains

    !> Initialize the atmosphere state with default values.
    subroutine atmosphere_state_init(state)
        type(atmosphere_state_t), intent(inout) :: state
        call atmosphere_state_reset_all(state)
    end subroutine atmosphere_state_init

    !> Finalize the atmosphere state, cleanup.
    subroutine atmosphere_state_finalize(state)
        type(atmosphere_state_t), intent(inout) :: state
        call atmosphere_state_reset_all(state)
    end subroutine atmosphere_state_finalize

    !> Reset all atmosphere state fields to canonical defaults.
    !!
    !! This routine is used by both initialization and finalization to avoid
    !! hidden carry-over state between model runs and instances.
    !! 
    !! In the initialization, all fields are initialized to zero or false by
    !! default in the type definition. This is generally sufficient, but explicitly
    !! resetting all fields ensures that any future additions to the state are
    !! properly initialized without relying on the default initialization behavior.
    !!
    !! In the finilize case this is a maximum safety solution, but can be adapted later.
    !! In practice, there are no allocated resources in the current state
    !! implementation, but this pattern allows for future-proofing if such resources are added later.
    !! 
    !!
    !! @param[inout] state Atmosphere state container
    subroutine atmosphere_state_reset_all(state)
        type(atmosphere_state_t), intent(inout) :: state

        ! Meteo input configuration
        state%swetr = 0
        state%swdivide = 0
        state%swmetdetail = 0
        state%swmeteo = 0
        state%swrain = 0
        state%swetsine = 0
        state%swinter = 0
        state%swMetFilAll = 0
        state%nmetdetail = 0
        state%altw = 0.0d0
        state%angstroma = 0.0d0
        state%angstromb = 0.0d0

        ! Current meteorological values
        state%tav = 0.0d0
        state%tavd = 0.0d0
        state%tmn = 0.0d0
        state%tmx = 0.0d0
        state%tmnr = 0.0d0
        state%rad = 0.0d0
        state%rh = 0.0d0
        state%lat = 0.0d0
        state%alt = 0.0d0
        state%daylp = 0.0d0

        ! Precipitation
        state%grai = 0.0d0
        state%graidt = 0.0d0
        state%nraida = 0.0d0
        state%nraidt = 0.0d0
        state%finterception = 0.0d0

        ! Evapotranspiration
        state%peva = 0.0d0
        state%pevaday = 0.0d0
        state%ptra = 0.0d0
        state%ptraday = 0.0d0
        state%tra = 0.0d0
        state%reva = 0.0d0
        state%atmdem = 0.0d0
        state%es0 = 0.0d0
        state%et0 = 0.0d0
        state%ew0 = 0.0d0

        ! Interception
        state%aintcdt = 0.0d0
        state%sicact = 0.0d0
        state%siccapact = 0.0d0

        ! Cumulative and intermediate fluxes
        call atmosphere_state_reset_cumulative(state)
        call atmosphere_state_reset_intermediate(state)

        ! Meteo file reading state
        state%daymeteo = 0
        state%yearmeteo = 0
        state%daynrfirst = 0
        state%daynrlast = 0
        state%wrecord = 0
        state%rainrec = 0
        state%timjan1 = 0.0d0
        state%metperiod = 0.0d0
        state%dtEventRain = 0.0d0

        ! Temperature history
        state%atmin7 = 0.0d0

        ! Flags
        state%fletsine = .false.
        state%flmeteodt = .false.
        state%flmetdetail = .false.
        state%flrainintens = .false.
        state%flupdmetdet = .false.

        ! ETSine and CN method persistent state
        state%tsunrise = 0.0d0
        state%tsunset = 0.0d0
        state%nod10_cn = 0
        state%icn_atm = 0
        state%z10_cn = 0.0d0

        ! Additional evaporation parameters
        state%empreva = 0.0d0
        state%fprecnosnow = 0.0d0

        ! Paths and files
        state%metfil = ''
        state%rainfil = ''
        state%pathatm = ''
    end subroutine atmosphere_state_reset_all

    !> Reset cumulative flux accumulators (called when flzerocumu=.true.)
    subroutine atmosphere_state_reset_cumulative(state)
        type(atmosphere_state_t), intent(inout) :: state
        
        state%cgrai = 0.0d0
        state%cnrai = 0.0d0
        state%caintc = 0.0d0
        state%cevap = 0.0d0
        state%cpeva = 0.0d0
        state%cptra = 0.0d0
    end subroutine atmosphere_state_reset_cumulative

    !> Reset intermediate flux accumulators (called when flzerointr=.true.)
    subroutine atmosphere_state_reset_intermediate(state)
        type(atmosphere_state_t), intent(inout) :: state
        
        state%inrai = 0.0d0
        state%igrai = 0.0d0
        state%iintc = 0.0d0
        state%ievap = 0.0d0
        state%ipeva = 0.0d0
        state%iptra = 0.0d0
        state%ies0 = 0.0d0
        state%iet0 = 0.0d0
        state%iew0 = 0.0d0
    end subroutine atmosphere_state_reset_intermediate

end module atmosphere_state_mod