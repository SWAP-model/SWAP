! ==============================================================================
! SWAP State Synchronization Module
! ==============================================================================
! This module provides bidirectional synchronization between the legacy
! module-level variables (in variables.f90) and the new state types
! (in swap_state_mod.f90).
!
! Purpose:
!   - Enable gradual refactoring from global variables to explicit state
!   - Maintain backwards compatibility during transition
!   - Allow multi-instance execution by snapshotting/restoring state
!
! Usage:
!   call state_from_variables(state)  ! Copy globals -> state (snapshot)
!   call state_to_variables(state)    ! Copy state -> globals (restore)
!
! During refactoring:
!   1. Original code continues to use module variables from 'variables'
!   2. New code works with explicit state types
!   3. Sync procedures bridge the two approaches
!   4. Eventually, variables.f90 can be deprecated
!
! Author: SWAP Development Team
! Date: 2026-01-31
! ==============================================================================

module swap_state_sync
    use swap_state_mod
    use swap_log, only: log_debug, log_info, to_str
    use variables  ! The legacy module with all global variables
    implicit none
    private

    ! Main synchronization procedures
    public :: state_from_variables
    public :: state_to_variables
    
    ! Component-level sync (for selective updates)
    public :: time_state_from_variables
    public :: time_state_to_variables
    public :: soil_state_from_variables
    public :: soil_state_to_variables
    public :: atmosphere_state_from_variables
    public :: atmosphere_state_to_variables
    public :: crop_state_from_variables
    public :: crop_state_to_variables
    public :: irrigation_state_from_variables
    public :: irrigation_state_to_variables
    public :: drainage_state_from_variables
    public :: drainage_state_to_variables

contains

    ! ==========================================================================
    ! Master Synchronization: Variables -> State
    ! ==========================================================================
    subroutine state_from_variables(state)
        type(swap_state_t), intent(inout) :: state
        
        call log_debug('sync', 'Syncing all variables -> state')
        
        call time_state_from_variables(state%time)
        call soil_state_from_variables(state%soil, state%numnod, state%numlay)
        call atmosphere_state_from_variables(state%atm)
        call crop_state_from_variables(state%crop, state%ncrop)
        call irrigation_state_from_variables(state%irrig)
        call drainage_state_from_variables(state%drain, state%nrlevs, state%numnod)
        
        call log_info('sync', 'State synchronized from variables')
    end subroutine state_from_variables

    ! ==========================================================================
    ! Master Synchronization: State -> Variables
    ! ==========================================================================
    subroutine state_to_variables(state)
        type(swap_state_t), intent(in) :: state
        
        call log_debug('sync', 'Syncing state -> all variables')
        
        call time_state_to_variables(state%time)
        call soil_state_to_variables(state%soil, state%numnod, state%numlay)
        call atmosphere_state_to_variables(state%atm)
        call crop_state_to_variables(state%crop, state%ncrop)
        call irrigation_state_to_variables(state%irrig)
        call drainage_state_to_variables(state%drain, state%nrlevs)
        
        call log_info('sync', 'Variables synchronized from state')
    end subroutine state_to_variables

    ! ==========================================================================
    ! Time State Synchronization
    ! ==========================================================================
    subroutine time_state_from_variables(tstate)
        type(time_state_t), intent(inout) :: tstate
        
        ! Time stepping
        tstate%dt = dt
        tstate%dtmax = dtmax
        tstate%dtmin = dtmin
        tstate%t = t
        tstate%t1900 = t1900
        tstate%tcum = tcum
        tstate%tend = tend
        tstate%tstart = tstart
        
        ! Counters
        tstate%daycum = daycum
        tstate%daynr = daynr
        tstate%imonth = imonth
        tstate%iyear = iyear
        tstate%iyearm1 = iyearm1
        tstate%isteps = isteps
        
        ! Output control
        tstate%ioutdat = ioutdat
        tstate%ioutdatint = ioutdatint
        tstate%cntper = cntper
        tstate%nprintday = nprintday
        tstate%nprintcount = nprintcount
        tstate%period = period
        
        ! Flags
        tstate%fldayend = fldayend
        tstate%fldaystart = fldaystart
        tstate%flrunend = flrunend
        tstate%flyearstart = flyearstart
        tstate%floutput = floutput
        tstate%flbaloutput = flbaloutput
        tstate%flheader = flheader
        tstate%flzerocumu = flzerocumu
        tstate%flzerointr = flzerointr
        tstate%fldecdt = fldecdt
        tstate%fldecdtmin = fldecdtmin
        tstate%fldtmin = fldtmin
        tstate%fldtreduce = fldtreduce
        tstate%flprintdt = flprintdt
        tstate%flprintshort = flprintshort
        tstate%floutputshort = floutputshort
        tstate%flSwapShared = flSwapShared
        
        ! Stress accumulation
        tstate%iqredwet_day = iqredwet_day
        tstate%iqreddry_day = iqreddry_day
        tstate%iqredsol_day = iqredsol_day
        tstate%iqredfrs_day = iqredfrs_day
        tstate%iptra_day = iptra_day
        
        ! Character variables
        tstate%date = date
        tstate%outfil = outfil
        tstate%pathwork = pathwork
        tstate%project = project
        tstate%swpfile = swpfile
        
        call log_debug('sync', 'time_state_from_variables: day=' // to_str(daynr) // ', year=' // to_str(iyear))
    end subroutine time_state_from_variables

    subroutine time_state_to_variables(tstate)
        type(time_state_t), intent(in) :: tstate
        
        ! Time stepping
        dt = tstate%dt
        dtmax = tstate%dtmax
        dtmin = tstate%dtmin
        t = tstate%t
        t1900 = tstate%t1900
        tcum = tstate%tcum
        tend = tstate%tend
        tstart = tstate%tstart
        
        ! Counters
        daycum = tstate%daycum
        daynr = tstate%daynr
        imonth = tstate%imonth
        iyear = tstate%iyear
        iyearm1 = tstate%iyearm1
        isteps = tstate%isteps
        
        ! Output control
        ioutdat = tstate%ioutdat
        ioutdatint = tstate%ioutdatint
        cntper = tstate%cntper
        nprintday = tstate%nprintday
        nprintcount = tstate%nprintcount
        period = tstate%period
        
        ! Flags
        fldayend = tstate%fldayend
        fldaystart = tstate%fldaystart
        flrunend = tstate%flrunend
        flyearstart = tstate%flyearstart
        floutput = tstate%floutput
        flbaloutput = tstate%flbaloutput
        flheader = tstate%flheader
        flzerocumu = tstate%flzerocumu
        flzerointr = tstate%flzerointr
        fldecdt = tstate%fldecdt
        fldecdtmin = tstate%fldecdtmin
        fldtmin = tstate%fldtmin
        fldtreduce = tstate%fldtreduce
        flprintdt = tstate%flprintdt
        flprintshort = tstate%flprintshort
        floutputshort = tstate%floutputshort
        flSwapShared = tstate%flSwapShared
        
        ! Stress accumulation
        iqredwet_day = tstate%iqredwet_day
        iqreddry_day = tstate%iqreddry_day
        iqredsol_day = tstate%iqredsol_day
        iqredfrs_day = tstate%iqredfrs_day
        iptra_day = tstate%iptra_day
        
        ! Character variables
        date = tstate%date
        outfil = tstate%outfil
        pathwork = tstate%pathwork
        project = tstate%project
        swpfile = tstate%swpfile
        
        call log_debug('sync', 'time_state_to_variables: day=' // to_str(tstate%daynr) // ', year=' // to_str(tstate%iyear))
    end subroutine time_state_to_variables

    ! ==========================================================================
    ! Soil State Synchronization
    ! ==========================================================================
    subroutine soil_state_from_variables(sstate, n_nod, n_lay)
        type(soil_state_t), intent(inout) :: sstate
        integer, intent(in) :: n_nod, n_lay
        integer :: i
        
        call log_debug('soil_sync', 'Starting soil_state_from_variables')
        
        ! Primary state variables (node-based)
        if (allocated(sstate%h) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(sstate%h))
                sstate%h(i) = h(i)
                sstate%theta(i) = theta(i)
                sstate%k(i) = k(i)
            end do
        end if
        
        ! Groundwater
        sstate%gwl = gwl
        sstate%nodgwl = nodgwl
        sstate%pond = pond
        
        ! Cumulative fluxes
        sstate%cqbot = cqbot
        sstate%crunoff = crunoff
        
        ! Scalars
        sstate%numnod = numnod
        sstate%numlay = numlay
        sstate%nsublay = nsublay
        sstate%numbit = numbit
        sstate%msteps = msteps
        sstate%CritDevh1Cp = CritDevh1Cp
        sstate%CritDevh2Cp = CritDevh2Cp
        sstate%FlRunoff = FlRunoff
        
        call log_debug('sync', 'soil_state_from_variables: gwl=' // to_str(real(gwl,4)))
    end subroutine soil_state_from_variables

    subroutine soil_state_to_variables(sstate, n_nod, n_lay)
        type(soil_state_t), intent(in) :: sstate
        integer, intent(in) :: n_nod, n_lay
        integer :: i
        
        ! Primary state variables (node-based)
        if (allocated(sstate%h) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(sstate%h))
                h(i) = sstate%h(i)
                theta(i) = sstate%theta(i)
                k(i) = sstate%k(i)
            end do
        end if
        
        ! Groundwater
        gwl = sstate%gwl
        nodgwl = sstate%nodgwl
        pond = sstate%pond
        
        ! Cumulative fluxes
        cqbot = sstate%cqbot
        crunoff = sstate%crunoff
        
        ! Scalars
        numnod = sstate%numnod
        numlay = sstate%numlay
        nsublay = sstate%nsublay
        numbit = sstate%numbit
        msteps = sstate%msteps
        CritDevh1Cp = sstate%CritDevh1Cp
        CritDevh2Cp = sstate%CritDevh2Cp
        FlRunoff = sstate%FlRunoff
        
        call log_debug('sync', 'soil_state_to_variables: gwl=' // to_str(real(sstate%gwl,4)))
    end subroutine soil_state_to_variables

    ! ==========================================================================
    ! Atmosphere State Synchronization
    ! ==========================================================================
    subroutine atmosphere_state_from_variables(astate)
        type(atmosphere_state_t), intent(inout) :: astate
        
        ! Current meteorological values
        astate%tav = tav
        astate%tavd = tavd
        astate%tmn = tmn
        astate%tmx = tmx
        astate%tmnr = tmnr
        astate%rad = rad
        astate%rh = rh
        astate%lat = lat
        astate%alt = alt
        astate%daylp = daylp
        
        ! Precipitation
        astate%grai = grai
        astate%graidt = graidt
        astate%nraida = nraida
        astate%nraidt = nraidt
        astate%finterception = finterception
        
        ! Evapotranspiration
        astate%peva = peva
        astate%pevaday = pevaday
        astate%ptra = ptra
        astate%ptraday = ptraday
        astate%tra = tra
        astate%atmdem = atmdem
        
        ! Cumulative values
        astate%cgrai = cgrai
        astate%cnrai = cnrai
        astate%caintc = caintc
        astate%cevap = cevap
        astate%cpeva = cpeva
        astate%cptra = cptra
        
        ! Intermediate values
        astate%inrai = inrai
        astate%ievap = ievap
        astate%ipeva = ipeva
        astate%iptra = iptra
        
        ! Meteo file reading state
        astate%daymeteo = daymeteo
        astate%yearmeteo = yearmeteo
        astate%daynrfirst = daynrfirst
        astate%daynrlast = daynrlast
        astate%wrecord = wrecord
        astate%rainrec = rainrec
        astate%timjan1 = timjan1
        astate%metperiod = metperiod
        astate%dtEventRain = dtEventRain
        
        ! Temperature history
        astate%atmin7 = atmin7
        
        ! Flags
        astate%fletsine = fletsine
        astate%flmeteodt = flmeteodt
        astate%flmetdetail = flmetdetail
        astate%flrainintens = flrainintens
        astate%flupdmetdet = flupdmetdet
        
        ! Paths
        astate%metfil = metfil
        astate%rainfil = rainfil
        astate%pathatm = pathatm
        
        call log_debug('sync', 'atmosphere_state_from_variables: tav=' // to_str(real(tav,4)))
    end subroutine atmosphere_state_from_variables

    subroutine atmosphere_state_to_variables(astate)
        type(atmosphere_state_t), intent(in) :: astate
        
        ! Current meteorological values
        tav = astate%tav
        tavd = astate%tavd
        tmn = astate%tmn
        tmx = astate%tmx
        tmnr = astate%tmnr
        rad = astate%rad
        rh = astate%rh
        lat = astate%lat
        alt = astate%alt
        daylp = astate%daylp
        
        ! Precipitation
        grai = astate%grai
        graidt = astate%graidt
        nraida = astate%nraida
        nraidt = astate%nraidt
        finterception = astate%finterception
        
        ! Evapotranspiration
        peva = astate%peva
        pevaday = astate%pevaday
        ptra = astate%ptra
        ptraday = astate%ptraday
        tra = astate%tra
        atmdem = astate%atmdem
        
        ! Cumulative values
        cgrai = astate%cgrai
        cnrai = astate%cnrai
        caintc = astate%caintc
        cevap = astate%cevap
        cpeva = astate%cpeva
        cptra = astate%cptra
        
        ! Intermediate values
        inrai = astate%inrai
        ievap = astate%ievap
        ipeva = astate%ipeva
        iptra = astate%iptra
        
        ! Meteo file reading state
        daymeteo = astate%daymeteo
        yearmeteo = astate%yearmeteo
        daynrfirst = astate%daynrfirst
        daynrlast = astate%daynrlast
        wrecord = astate%wrecord
        rainrec = astate%rainrec
        timjan1 = astate%timjan1
        metperiod = astate%metperiod
        dtEventRain = astate%dtEventRain
        
        ! Temperature history
        atmin7 = astate%atmin7
        
        ! Flags
        fletsine = astate%fletsine
        flmeteodt = astate%flmeteodt
        flmetdetail = astate%flmetdetail
        flrainintens = astate%flrainintens
        flupdmetdet = astate%flupdmetdet
        
        ! Paths
        metfil = astate%metfil
        rainfil = astate%rainfil
        pathatm = astate%pathatm
        
        call log_debug('sync', 'atmosphere_state_to_variables: tav=' // to_str(real(astate%tav,4)))
    end subroutine atmosphere_state_to_variables

    ! ==========================================================================
    ! Crop State Synchronization
    ! ==========================================================================
    subroutine crop_state_from_variables(cstate, n_crop)
        type(crop_state_t), intent(inout) :: cstate
        integer, intent(in) :: n_crop
        
        ! Development
        cstate%dvs = dvs
        cstate%dvsend = dvsend
        cstate%tsum = tsum
        cstate%tsumgerm = tsumgerm
        cstate%daycrop = daycrop
        cstate%daygrowth = daygrowth
        cstate%icrop = icrop
        
        ! Leaf area and canopy
        cstate%lai = lai
        cstate%laipot = laipot
        cstate%laimax = laimax
        cstate%laiexp = laiexp
        cstate%gc = gc
        cstate%ch = ch
        cstate%cf = cf
        cstate%albedo = albedo
        
        ! Biomass (actual)
        cstate%wlv = wlv
        cstate%wst = wst
        cstate%wrt = wrt
        cstate%wso = wso
        cstate%cwdm = cwdm
        cstate%tadw = tadw
        
        ! Biomass (potential)
        cstate%wlvpot = wlvpot
        cstate%wstpot = wstpot
        cstate%wrtpot = wrtpot
        cstate%wsopot = wsopot
        cstate%cwdmpot = cwdmpot
        cstate%tadwpot = tadwpot
        
        ! Dead matter
        cstate%dwlv = dwlv
        cstate%dwst = dwst
        cstate%dwrt = dwrt
        cstate%dwso = dwso
        
        ! Assimilation
        cstate%gasst = gasst
        cstate%gasstpot = gasstpot
        cstate%mrest = mrest
        cstate%mrestpot = mrestpot
        
        ! Rooting
        cstate%rd = rd
        cstate%rdpot = rdpot
        cstate%rdi = rdi
        cstate%rri = rri
        cstate%rdc = rdc
        cstate%rdm = rdm
        cstate%rdmax = rdmax
        cstate%noddrz = noddrz
        
        ! Root uptake
        cstate%hlim1 = hlim1
        cstate%hlim2u = hlim2u
        cstate%hlim2l = hlim2l
        cstate%hlim3h = hlim3h
        cstate%hlim3l = hlim3l
        cstate%hlim4 = hlim4
        
        ! Stress
        cstate%reltr = reltr
        cstate%alphacrit = alphacrit
        cstate%alpJvLier = alpJvLier
        
        ! Grass
        cstate%tagp = tagp
        cstate%tagppot = tagppot
        cstate%tagpt = tagpt
        cstate%mowrest = mowrest
        cstate%cuptgraz = cuptgraz
        cstate%iharvest = iharvest
        cstate%iseqgm = iseqgm
        
        ! Flags
        cstate%flCropCalendar = flCropCalendar
        cstate%flCropEmergence = flCropEmergence
        cstate%flCropHarvest = flCropHarvest
        cstate%flanthesis = flanthesis
        cstate%flHarvest = flHarvest
        cstate%flHarvestDay = flHarvestDay
        cstate%flGrazing = flGrazing
        cstate%flCO2 = flCO2
        
        ! CO2 factors
        cstate%fco2amax = fco2amax
        cstate%fco2eff = fco2eff
        cstate%fco2tra = fco2tra
        
        ! Path
        cstate%pathcrop = pathcrop
        
        call log_debug('sync', 'crop_state_from_variables: dvs=' // to_str(real(dvs,4)) // ', lai=' // to_str(real(lai,4)))
    end subroutine crop_state_from_variables

    subroutine crop_state_to_variables(cstate, n_crop)
        type(crop_state_t), intent(in) :: cstate
        integer, intent(in) :: n_crop
        
        ! Development
        dvs = cstate%dvs
        dvsend = cstate%dvsend
        tsum = cstate%tsum
        tsumgerm = cstate%tsumgerm
        daycrop = cstate%daycrop
        daygrowth = cstate%daygrowth
        icrop = cstate%icrop
        
        ! Leaf area and canopy
        lai = cstate%lai
        laipot = cstate%laipot
        laimax = cstate%laimax
        laiexp = cstate%laiexp
        gc = cstate%gc
        ch = cstate%ch
        cf = cstate%cf
        albedo = cstate%albedo
        
        ! Biomass (actual)
        wlv = cstate%wlv
        wst = cstate%wst
        wrt = cstate%wrt
        wso = cstate%wso
        cwdm = cstate%cwdm
        tadw = cstate%tadw
        
        ! Biomass (potential)
        wlvpot = cstate%wlvpot
        wstpot = cstate%wstpot
        wrtpot = cstate%wrtpot
        wsopot = cstate%wsopot
        cwdmpot = cstate%cwdmpot
        tadwpot = cstate%tadwpot
        
        ! Dead matter
        dwlv = cstate%dwlv
        dwst = cstate%dwst
        dwrt = cstate%dwrt
        dwso = cstate%dwso
        
        ! Assimilation
        gasst = cstate%gasst
        gasstpot = cstate%gasstpot
        mrest = cstate%mrest
        mrestpot = cstate%mrestpot
        
        ! Rooting
        rd = cstate%rd
        rdpot = cstate%rdpot
        rdi = cstate%rdi
        rri = cstate%rri
        rdc = cstate%rdc
        rdm = cstate%rdm
        rdmax = cstate%rdmax
        noddrz = cstate%noddrz
        
        ! Root uptake
        hlim1 = cstate%hlim1
        hlim2u = cstate%hlim2u
        hlim2l = cstate%hlim2l
        hlim3h = cstate%hlim3h
        hlim3l = cstate%hlim3l
        hlim4 = cstate%hlim4
        
        ! Stress
        reltr = cstate%reltr
        alphacrit = cstate%alphacrit
        alpJvLier = cstate%alpJvLier
        
        ! Grass
        tagp = cstate%tagp
        tagppot = cstate%tagppot
        tagpt = cstate%tagpt
        mowrest = cstate%mowrest
        cuptgraz = cstate%cuptgraz
        iharvest = cstate%iharvest
        iseqgm = cstate%iseqgm
        
        ! Flags
        flCropCalendar = cstate%flCropCalendar
        flCropEmergence = cstate%flCropEmergence
        flCropHarvest = cstate%flCropHarvest
        flanthesis = cstate%flanthesis
        flHarvest = cstate%flHarvest
        flHarvestDay = cstate%flHarvestDay
        flGrazing = cstate%flGrazing
        flCO2 = cstate%flCO2
        
        ! CO2 factors
        fco2amax = cstate%fco2amax
        fco2eff = cstate%fco2eff
        fco2tra = cstate%fco2tra
        
        ! Path
        pathcrop = cstate%pathcrop
        
        call log_debug('sync', 'crop_state_to_variables: dvs=' // to_str(real(cstate%dvs,4)) // ', lai=' // to_str(real(cstate%lai,4)))
    end subroutine crop_state_to_variables

    ! ==========================================================================
    ! Irrigation State Synchronization
    ! ==========================================================================
    subroutine irrigation_state_from_variables(istate)
        type(irrigation_state_t), intent(inout) :: istate
        
        istate%gird = gird
        istate%nird = nird
        istate%cirrs = cirrs
        istate%igird = igird
        istate%inird = inird
        
        istate%irrigevent = irrigevent
        istate%nirri = nirri
        istate%isua = isua
        istate%dayfix = dayfix
        
        istate%flirrigate = flirrigate
        istate%flheadirg = flheadirg
        
        call log_debug('sync', 'irrigation_state_from_variables')
    end subroutine irrigation_state_from_variables

    subroutine irrigation_state_to_variables(istate)
        type(irrigation_state_t), intent(in) :: istate
        
        gird = istate%gird
        nird = istate%nird
        cirrs = istate%cirrs
        igird = istate%igird
        inird = istate%inird
        
        irrigevent = istate%irrigevent
        nirri = istate%nirri
        isua = istate%isua
        dayfix = istate%dayfix
        
        flirrigate = istate%flirrigate
        flheadirg = istate%flheadirg
        
        call log_debug('sync', 'irrigation_state_to_variables')
    end subroutine irrigation_state_to_variables

    ! ==========================================================================
    ! Drainage State Synchronization
    ! ==========================================================================
    subroutine drainage_state_from_variables(dstate, n_levs, n_nod)
        type(drainage_state_t), intent(inout) :: dstate
        integer, intent(in) :: n_levs, n_nod
        integer :: i, j
        
        ! Copy drainage fluxes per level
        if (allocated(dstate%qdrain) .and. n_levs > 0) then
            do i = 1, min(n_levs, size(dstate%qdrain))
                dstate%qdrain(i) = qdrain(i)
            end do
        end if
        
        ! Scalars
        dstate%qdrtot = qdrtot
        dstate%nrlevs = nrlevs
        dstate%dramet = dramet
        dstate%fldrain = fldrain
        
        call log_debug('sync', 'drainage_state_from_variables: nrlevs=' // to_str(nrlevs))
    end subroutine drainage_state_from_variables

    subroutine drainage_state_to_variables(dstate, n_levs)
        type(drainage_state_t), intent(in) :: dstate
        integer, intent(in) :: n_levs
        integer :: i
        
        ! Copy drainage fluxes per level
        if (allocated(dstate%qdrain) .and. n_levs > 0) then
            do i = 1, min(n_levs, size(dstate%qdrain))
                qdrain(i) = dstate%qdrain(i)
            end do
        end if
        
        ! Scalars
        qdrtot = dstate%qdrtot
        nrlevs = dstate%nrlevs
        dramet = dstate%dramet
        fldrain = dstate%fldrain
        
        call log_debug('sync', 'drainage_state_to_variables: nrlevs=' // to_str(dstate%nrlevs))
    end subroutine drainage_state_to_variables

end module swap_state_sync
