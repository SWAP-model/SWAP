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
    public :: surfacewater_state_from_variables
    public :: surfacewater_state_to_variables
    public :: boundary_state_from_variables
    public :: boundary_state_to_variables
    public :: solute_state_from_variables
    public :: solute_state_to_variables
    public :: heat_state_from_variables
    public :: heat_state_to_variables

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
        call boundary_state_from_variables(state%boundary, MADAY)
        call solute_state_from_variables(state%solute, state%numnod, state%numlay, state%nrlevs)
        call heat_state_from_variables(state%heat, state%numnod, state%numlay)
        
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
        call boundary_state_to_variables(state%boundary, MADAY)
        call solute_state_to_variables(state%solute, state%numnod, state%numlay, state%nrlevs)
        call heat_state_to_variables(state%heat, state%numnod, state%numlay)
        
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
            call log_debug('soil_sync', 'Copied h, theta, k arrays for ' // to_str(n_nod) // ' nodes')
        end if
        
        ! Previous timestep values
        if (allocated(sstate%hm1) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(sstate%hm1))
                sstate%hm1(i) = hm1(i)
                sstate%thetm1(i) = thetm1(i)
            end do
        end if
        
        ! Fluxes
        if (allocated(sstate%q) .and. n_nod > 0) then
            do i = 1, min(n_nod+1, size(sstate%q))
                sstate%q(i) = q(i)
            end do
        end if
        
        if (allocated(sstate%qrot)) then
            do i = 1, min(n_nod, size(sstate%qrot))
                sstate%qrot(i) = qrot(i)
            end do
        end if
        
        ! Groundwater
        sstate%gwl = gwl
        sstate%gwlm1 = gwlm1
        sstate%gwli = gwli
        sstate%gwlinp = gwlinp
        sstate%nodgwl = nodgwl
        sstate%npegwl = npegwl
        sstate%bpegwl = bpegwl
        sstate%fllowgwl = fllowgwl
        
        ! Surface/ponding
        sstate%pond = pond
        sstate%pondm1 = pondm1
        sstate%pondmx = pondmx
        sstate%qtop = qtop
        sstate%qbot = qbot
        
        ! Cumulative fluxes
        sstate%cqbot = cqbot
        sstate%cqbotdo = cqbotdo
        sstate%cqbotup = cqbotup
        sstate%cqtdo = cqtdo
        sstate%cqtup = cqtup
        sstate%cqrot = cqrot
        sstate%cqdra = cqdra
        sstate%crunoff = crunoff
        sstate%crunon = crunon
        
        ! Storage
        sstate%volact = volact
        sstate%volini = volini
        sstate%volm1 = volm1
        
        ! Evaporation reduction (Boesten/Black)
        sstate%saev = saev
        sstate%spev = spev
        sstate%ldwet = ldwet
        sstate%cofred = cofred
        
        ! Iteration control
        sstate%numnod = numnod
        sstate%numlay = numlay
        sstate%nsublay = nsublay
        sstate%numbit = numbit
        sstate%msteps = msteps
        sstate%CritDevh1Cp = CritDevh1Cp
        sstate%CritDevh2Cp = CritDevh2Cp
        sstate%CritDevMasBal = CritDevMasBal
        sstate%gwlconv = gwlconv
        
        ! Flags
        sstate%FlRunoff = FlRunoff
        sstate%fldrain = fldrain
        sstate%swhyst = swhyst
        
        ! Headcalc iteration tracking (formerly SAVE variables)
        sstate%flwarn_hc = flwarn_hc
        sstate%iwarn_hc = iwarn_hc
        sstate%nstep_hc = nstep_hc
        
        call log_debug('soil_sync', 'soil_state_from_variables complete: gwl=' // &
                       to_str(real(gwl,4)) // ', pond=' // to_str(real(pond,4)))
    end subroutine soil_state_from_variables

    subroutine soil_state_to_variables(sstate, n_nod, n_lay)
        type(soil_state_t), intent(in) :: sstate
        integer, intent(in) :: n_nod, n_lay
        integer :: i
        
        call log_debug('soil_sync', 'Starting soil_state_to_variables')
        
        ! Primary state variables (node-based)
        if (allocated(sstate%h) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(sstate%h))
                h(i) = sstate%h(i)
                theta(i) = sstate%theta(i)
                k(i) = sstate%k(i)
            end do
            call log_debug('soil_sync', 'Restored h, theta, k arrays for ' // to_str(n_nod) // ' nodes')
        end if
        
        ! Previous timestep values
        if (allocated(sstate%hm1) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(sstate%hm1))
                hm1(i) = sstate%hm1(i)
                thetm1(i) = sstate%thetm1(i)
            end do
        end if
        
        ! Fluxes
        if (allocated(sstate%q) .and. n_nod > 0) then
            do i = 1, min(n_nod+1, size(sstate%q))
                q(i) = sstate%q(i)
            end do
        end if
        
        if (allocated(sstate%qrot)) then
            do i = 1, min(n_nod, size(sstate%qrot))
                qrot(i) = sstate%qrot(i)
            end do
        end if
        
        ! Groundwater
        gwl = sstate%gwl
        gwlm1 = sstate%gwlm1
        gwli = sstate%gwli
        gwlinp = sstate%gwlinp
        nodgwl = sstate%nodgwl
        npegwl = sstate%npegwl
        bpegwl = sstate%bpegwl
        fllowgwl = sstate%fllowgwl
        
        ! Surface/ponding
        pond = sstate%pond
        pondm1 = sstate%pondm1
        pondmx = sstate%pondmx
        qtop = sstate%qtop
        qbot = sstate%qbot
        
        ! Cumulative fluxes
        cqbot = sstate%cqbot
        cqbotdo = sstate%cqbotdo
        cqbotup = sstate%cqbotup
        cqtdo = sstate%cqtdo
        cqtup = sstate%cqtup
        cqrot = sstate%cqrot
        cqdra = sstate%cqdra
        crunoff = sstate%crunoff
        crunon = sstate%crunon
        
        ! Storage
        volact = sstate%volact
        volini = sstate%volini
        volm1 = sstate%volm1
        
        ! Evaporation reduction (Boesten/Black)
        saev = sstate%saev
        spev = sstate%spev
        ldwet = sstate%ldwet
        cofred = sstate%cofred
        
        ! Iteration control
        numnod = sstate%numnod
        numlay = sstate%numlay
        nsublay = sstate%nsublay
        numbit = sstate%numbit
        msteps = sstate%msteps
        CritDevh1Cp = sstate%CritDevh1Cp
        CritDevh2Cp = sstate%CritDevh2Cp
        CritDevMasBal = sstate%CritDevMasBal
        gwlconv = sstate%gwlconv
        
        ! Flags
        FlRunoff = sstate%FlRunoff
        fldrain = sstate%fldrain
        swhyst = sstate%swhyst
        
        ! Headcalc iteration tracking (formerly SAVE variables)
        flwarn_hc = sstate%flwarn_hc
        iwarn_hc = sstate%iwarn_hc
        nstep_hc = sstate%nstep_hc
        
        call log_debug('soil_sync', 'soil_state_to_variables complete: gwl=' // &
                       to_str(real(sstate%gwl,4)) // ', pond=' // to_str(real(sstate%pond,4)))
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
        
        ! ! ETSine SAVE variables (from meteodt.f90)
        ! astate%tsunrise = tsunrise_atm
        ! astate%tsunset = tsunset_atm
        
        ! ! CN runoff method SAVE variables (from meteoday.f90)
        ! astate%nod10_cn = nod10_cn
        ! astate%icn = icn_atm
        ! astate%z10_cn = z10_cn
        
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
        
        ! ! ETSine SAVE variables (from meteodt.f90)
        ! tsunrise_atm = astate%tsunrise
        ! tsunset_atm = astate%tsunset
        
        ! ! CN runoff method SAVE variables (from meteoday.f90)
        ! nod10_cn = astate%nod10_cn
        ! icn_atm = astate%icn
        ! z10_cn = astate%z10_cn
        
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
                dstate%cqdrain(i) = cqdrain(i)
                dstate%cqdrainin(i) = cqdrainin(i)
                dstate%cqdrainout(i) = cqdrainout(i)
                dstate%drainl(i) = drainl(i)
                dstate%drares(i) = drares(i)
                dstate%infres(i) = infres(i)
                dstate%L(i) = l(i)
                dstate%wetper(i) = wetper(i)
                dstate%zbotdr(i) = zbotdr(i)
                dstate%rdrain(i) = rdrain(i)
                dstate%rinfi(i) = rinfi(i)
                dstate%rentry(i) = rentry(i)
                dstate%rexit(i) = rexit(i)
                dstate%gwlinf(i) = gwlinf(i)
                dstate%widthr(i) = widthr(i)
                dstate%taludr(i) = taludr(i)
                dstate%swallo(i) = swallo(i)
                dstate%swdtyp(i) = swdtyp(i)
                dstate%swtopdislay(i) = swtopdislay(i)
                dstate%zTopDisLay(i) = zTopDisLay(i)
                dstate%fTopDisLay(i) = fTopDisLay(i)
            end do
        end if
        
        ! Spatial arrays
        if (allocated(dstate%qdra) .and. n_levs > 0 .and. n_nod > 0) then
            do i = 1, min(n_levs, size(dstate%qdra, 1))
                do j = 1, min(n_nod, size(dstate%qdra, 2))
                    dstate%qdra(i,j) = qdra(i,j)
                    dstate%inqdra(i,j) = inqdra(i,j)
                    dstate%inqdra_in(i,j) = inqdra_in(i,j)
                    dstate%inqdra_out(i,j) = inqdra_out(i,j)
                end do
            end do
        end if
        
        if (allocated(dstate%qdraincomp) .and. n_nod > 0) then
            do j = 1, min(n_nod, size(dstate%qdraincomp))
                dstate%qdraincomp(j) = qdraincomp(j)
            end do
        end if
        
        ! Scalars
        dstate%qdrtot = qdrtot
        dstate%iqdra = iqdra
        dstate%cqdra = cqdra
        dstate%nrlevs = nrlevs
        dstate%nrpri = nrpri
        dstate%dramet = dramet
        dstate%swdivd = swdivd
        dstate%swdislay = swdislay
        dstate%basegw = basegw
        dstate%entres = entres
        dstate%shape = shape
        
        ! Interflow parameters
        dstate%cofintfl = cofintfl
        dstate%expintfl = expintfl
        dstate%swnrsrf = swnrsrf
        dstate%SwTopnrsrf = SwTopnrsrf
        dstate%rsurfdeep = rsurfdeep
        dstate%rsurfshallow = rsurfshallow
        dstate%FacDpthInf = FacDpthInf
        dstate%Swdivdinf = Swdivdinf
        
        dstate%fldrain = fldrain
        
        call log_debug('sync', 'drainage_state_from_variables: nrlevs=' // to_str(nrlevs))
    end subroutine drainage_state_from_variables

    subroutine drainage_state_to_variables(dstate, n_levs)
        type(drainage_state_t), intent(in) :: dstate
        integer, intent(in) :: n_levs
        integer :: i, j, n_nod
        
        ! Determine number of nodes from array size
        n_nod = 0
        if (allocated(dstate%qdra)) n_nod = size(dstate%qdra, 2)
        
        ! Copy drainage fluxes per level
        if (allocated(dstate%qdrain) .and. n_levs > 0) then
            do i = 1, min(n_levs, size(dstate%qdrain))
                qdrain(i) = dstate%qdrain(i)
                cqdrain(i) = dstate%cqdrain(i)
                cqdrainin(i) = dstate%cqdrainin(i)
                cqdrainout(i) = dstate%cqdrainout(i)
                drainl(i) = dstate%drainl(i)
                drares(i) = dstate%drares(i)
                infres(i) = dstate%infres(i)
                l(i) = dstate%L(i)
                wetper(i) = dstate%wetper(i)
                zbotdr(i) = dstate%zbotdr(i)
                rdrain(i) = dstate%rdrain(i)
                rinfi(i) = dstate%rinfi(i)
                rentry(i) = dstate%rentry(i)
                rexit(i) = dstate%rexit(i)
                gwlinf(i) = dstate%gwlinf(i)
                widthr(i) = dstate%widthr(i)
                taludr(i) = dstate%taludr(i)
                swallo(i) = dstate%swallo(i)
                swdtyp(i) = dstate%swdtyp(i)
                swtopdislay(i) = dstate%swtopdislay(i)
                zTopDisLay(i) = dstate%zTopDisLay(i)
                fTopDisLay(i) = dstate%fTopDisLay(i)
            end do
        end if
        
        ! Spatial arrays
        if (allocated(dstate%qdra) .and. n_levs > 0 .and. n_nod > 0) then
            do i = 1, min(n_levs, size(dstate%qdra, 1))
                do j = 1, min(n_nod, size(dstate%qdra, 2))
                    qdra(i,j) = dstate%qdra(i,j)
                    inqdra(i,j) = dstate%inqdra(i,j)
                    inqdra_in(i,j) = dstate%inqdra_in(i,j)
                    inqdra_out(i,j) = dstate%inqdra_out(i,j)
                end do
            end do
        end if
        
        if (allocated(dstate%qdraincomp) .and. n_nod > 0) then
            do j = 1, min(n_nod, size(dstate%qdraincomp))
                qdraincomp(j) = dstate%qdraincomp(j)
            end do
        end if
        
        ! Scalars
        qdrtot = dstate%qdrtot
        iqdra = dstate%iqdra
        cqdra = dstate%cqdra
        nrlevs = dstate%nrlevs
        nrpri = dstate%nrpri
        dramet = dstate%dramet
        swdivd = dstate%swdivd
        swdislay = dstate%swdislay
        basegw = dstate%basegw
        entres = dstate%entres
        shape = dstate%shape
        
        ! Interflow parameters
        cofintfl = dstate%cofintfl
        expintfl = dstate%expintfl
        swnrsrf = dstate%swnrsrf
        SwTopnrsrf = dstate%SwTopnrsrf
        rsurfdeep = dstate%rsurfdeep
        rsurfshallow = dstate%rsurfshallow
        FacDpthInf = dstate%FacDpthInf
        Swdivdinf = dstate%Swdivdinf
        
        fldrain = dstate%fldrain
        
        call log_debug('sync', 'drainage_state_to_variables: nrlevs=' // to_str(dstate%nrlevs))
    end subroutine drainage_state_to_variables

    ! ==========================================================================
    ! Surface Water State Synchronization
    ! ==========================================================================
    subroutine surfacewater_state_from_variables(swstate, n_mper, n_levs)
        type(surfacewater_state_t), intent(inout) :: swstate
        integer, intent(in) :: n_mper, n_levs
        integer :: i, j
        
        ! Water levels
        swstate%wlp = wlp
        swstate%wls = wls
        swstate%wlsold = wlsold
        swstate%wlstar = wlstar
        swstate%hwlman = hwlman
        swstate%vtair = vtair
        
        ! Water level history
        if (allocated(swstate%wlsbak)) then
            do i = 1, min(4, size(swstate%wlsbak))
                swstate%wlsbak(i) = wlsbak(i)
            end do
        end if
        
        ! Storage and fluxes
        swstate%swst = swst
        swstate%swstini = swstini
        swstate%qdrd = qdrd
        swstate%cqdrd = cqdrd
        swstate%cwsupp = cwsupp
        swstate%cwout = cwout
        swstate%runots = runots
        swstate%QRapDra = QRapDra
        
        ! Management
        swstate%imper = imper
        swstate%nmper = nmper
        swstate%numadj = numadj
        swstate%swsrf = swsrf
        swstate%swsec = swsec
        swstate%swqhr = swqhr
        swstate%osswlm = osswlm
        swstate%overfl = overfl
        swstate%fldecdt = fldecdt
        swstate%fldtmin = fldtmin
        
        call log_debug('sync', 'surfacewater_state_from_variables: nmper=' // to_str(nmper))
    end subroutine surfacewater_state_from_variables

    subroutine surfacewater_state_to_variables(swstate, n_mper)
        type(surfacewater_state_t), intent(in) :: swstate
        integer, intent(in) :: n_mper
        integer :: i
        
        ! Water levels
        wlp = swstate%wlp
        wls = swstate%wls
        wlsold = swstate%wlsold
        wlstar = swstate%wlstar
        hwlman = swstate%hwlman
        vtair = swstate%vtair
        
        ! Water level history
        if (allocated(swstate%wlsbak)) then
            do i = 1, min(4, size(swstate%wlsbak))
                wlsbak(i) = swstate%wlsbak(i)
            end do
        end if
        
        ! Storage and fluxes
        swst = swstate%swst
        swstini = swstate%swstini
        qdrd = swstate%qdrd
        cqdrd = swstate%cqdrd
        cwsupp = swstate%cwsupp
        cwout = swstate%cwout
        runots = swstate%runots
        QRapDra = swstate%QRapDra
        
        ! Management
        imper = swstate%imper
        nmper = swstate%nmper
        numadj = swstate%numadj
        swsrf = swstate%swsrf
        swsec = swstate%swsec
        swqhr = swstate%swqhr
        osswlm = swstate%osswlm
        overfl = swstate%overfl
        fldecdt = swstate%fldecdt
        fldtmin = swstate%fldtmin
        
        call log_debug('sync', 'surfacewater_state_to_variables: nmper=' // to_str(swstate%nmper))
    end subroutine surfacewater_state_to_variables

    ! ==========================================================================
    ! Boundary State Synchronization
    ! ==========================================================================
    subroutine boundary_state_from_variables(bstate, n_day)
        type(boundary_state_t), intent(inout) :: bstate
        integer, intent(in) :: n_day
        integer :: i
        
        ! Bottom boundary - configuration
        bstate%swbotb = swbotb
        bstate%swbotb3Impl = swbotb3Impl
        bstate%SwBotb3ResVert = SwBotb3ResVert
        bstate%swqhbot = swqhbot
        bstate%swcofqhc = swcofqhc
        bstate%sw2 = sw2
        bstate%sw3 = sw3
        bstate%sw4 = sw4
        
        ! Bottom boundary - values
        bstate%qbot = qbot
        bstate%qbot_nonfrozen = qbot_nonfrozen
        bstate%hbot = hbot
        bstate%iqbot = iqbot
        bstate%cqbot = cqbot
        bstate%cqbotdo = cqbotdo
        bstate%cqbotup = cqbotup
        bstate%deepgw = deepgw
        
        ! Aquifer parameters
        bstate%aqave = aqave
        bstate%aqamp = aqamp
        bstate%aqper = aqper
        bstate%aqtmax = aqtmax
        bstate%rimlay = rimlay
        bstate%hdrain = hdrain
        bstate%shape = shape
        
        ! Sine function parameters
        bstate%sinave = sinave
        bstate%sinamp = sinamp
        bstate%sinmax = sinmax
        
        ! Flux-head relationships
        bstate%cofqha = cofqha
        bstate%cofqhb = cofqhb
        bstate%cofqhc = cofqhc
        
        ! Lysimeter
        bstate%hplate = hplate
        
        ! Tables - copy if allocated
        if (allocated(bstate%gwltab)) then
            do i = 1, min(size(bstate%gwltab), size(gwltab))
                bstate%gwltab(i) = gwltab(i)
            end do
        end if
        if (allocated(bstate%haqtab)) then
            do i = 1, min(size(bstate%haqtab), size(haqtab))
                bstate%haqtab(i) = haqtab(i)
            end do
        end if
        if (allocated(bstate%qbotab)) then
            do i = 1, min(size(bstate%qbotab), size(qbotab))
                bstate%qbotab(i) = qbotab(i)
            end do
        end if
        if (allocated(bstate%hbotab)) then
            do i = 1, min(size(bstate%hbotab), size(hbotab))
                bstate%hbotab(i) = hbotab(i)
            end do
        end if
        
        ! Top boundary - configuration
        bstate%swpondmx = swpondmx
        bstate%swredu = swredu
        
        ! Top boundary - values
        bstate%pondmx = pondmx
        bstate%hatm = hatm
        bstate%hsurf = hsurf
        bstate%rsro = rsro
        bstate%rsroexp = rsroexp
        bstate%runon = runon
        bstate%runots = runots
        bstate%crunoff = crunoff
        bstate%crunon = crunon
        bstate%iruno = iruno
        bstate%irunon = irunon
        
        ! Surface/ponding
        bstate%qtop = qtop
        bstate%q0 = q0
        bstate%h0max = h0max
        bstate%k1max = k1max
        bstate%QMpLatSs = QMpLatSs
        
        ! Runon array
        if (allocated(bstate%runonarr)) then
            do i = 1, min(size(bstate%runonarr), n_day)
                bstate%runonarr(i) = runonarr(i)
            end do
        end if
        
        ! Ponding table
        if (allocated(bstate%pondmxtab)) then
            do i = 1, min(size(bstate%pondmxtab), size(pondmxtab))
                bstate%pondmxtab(i) = pondmxtab(i)
            end do
        end if
        
        ! Inundation
        bstate%cinund = cinund
        
        ! Flags
        bstate%FlRunoff = FlRunoff
        bstate%flrunon = flrunon
        bstate%ftoph = ftoph
        
        call log_debug('sync', 'boundary_state_from_variables: swbotb=' // to_str(swbotb))
    end subroutine boundary_state_from_variables

    subroutine boundary_state_to_variables(bstate, n_day)
        type(boundary_state_t), intent(in) :: bstate
        integer, intent(in) :: n_day
        integer :: i
        
        ! Bottom boundary - configuration
        swbotb = bstate%swbotb
        swbotb3Impl = bstate%swbotb3Impl
        SwBotb3ResVert = bstate%SwBotb3ResVert
        swqhbot = bstate%swqhbot
        swcofqhc = bstate%swcofqhc
        sw2 = bstate%sw2
        sw3 = bstate%sw3
        sw4 = bstate%sw4
        
        ! Bottom boundary - values
        qbot = bstate%qbot
        qbot_nonfrozen = bstate%qbot_nonfrozen
        hbot = bstate%hbot
        iqbot = bstate%iqbot
        cqbot = bstate%cqbot
        cqbotdo = bstate%cqbotdo
        cqbotup = bstate%cqbotup
        deepgw = bstate%deepgw
        
        ! Aquifer parameters
        aqave = bstate%aqave
        aqamp = bstate%aqamp
        aqper = bstate%aqper
        aqtmax = bstate%aqtmax
        rimlay = bstate%rimlay
        hdrain = bstate%hdrain
        shape = bstate%shape
        
        ! Sine function parameters
        sinave = bstate%sinave
        sinamp = bstate%sinamp
        sinmax = bstate%sinmax
        
        ! Flux-head relationships
        cofqha = bstate%cofqha
        cofqhb = bstate%cofqhb
        cofqhc = bstate%cofqhc
        
        ! Lysimeter
        hplate = bstate%hplate
        
        ! Tables - copy if allocated
        if (allocated(bstate%gwltab)) then
            do i = 1, min(size(bstate%gwltab), size(gwltab))
                gwltab(i) = bstate%gwltab(i)
            end do
        end if
        if (allocated(bstate%haqtab)) then
            do i = 1, min(size(bstate%haqtab), size(haqtab))
                haqtab(i) = bstate%haqtab(i)
            end do
        end if
        if (allocated(bstate%qbotab)) then
            do i = 1, min(size(bstate%qbotab), size(qbotab))
                qbotab(i) = bstate%qbotab(i)
            end do
        end if
        if (allocated(bstate%hbotab)) then
            do i = 1, min(size(bstate%hbotab), size(hbotab))
                hbotab(i) = bstate%hbotab(i)
            end do
        end if
        
        ! Top boundary - configuration
        swpondmx = bstate%swpondmx
        swredu = bstate%swredu
        
        ! Top boundary - values
        pondmx = bstate%pondmx
        hatm = bstate%hatm
        hsurf = bstate%hsurf
        rsro = bstate%rsro
        rsroexp = bstate%rsroexp
        runon = bstate%runon
        runots = bstate%runots
        crunoff = bstate%crunoff
        crunon = bstate%crunon
        iruno = bstate%iruno
        irunon = bstate%irunon
        
        ! Surface/ponding
        qtop = bstate%qtop
        q0 = bstate%q0
        h0max = bstate%h0max
        k1max = bstate%k1max
        QMpLatSs = bstate%QMpLatSs
        
        ! Runon array
        if (allocated(bstate%runonarr)) then
            do i = 1, min(size(bstate%runonarr), n_day)
                runonarr(i) = bstate%runonarr(i)
            end do
        end if
        
        ! Ponding table
        if (allocated(bstate%pondmxtab)) then
            do i = 1, min(size(bstate%pondmxtab), size(pondmxtab))
                pondmxtab(i) = bstate%pondmxtab(i)
            end do
        end if
        
        ! Inundation
        cinund = bstate%cinund
        
        ! Flags
        FlRunoff = bstate%FlRunoff
        flrunon = bstate%flrunon
        ftoph = bstate%ftoph
        
        call log_debug('sync', 'boundary_state_to_variables: swbotb=' // to_str(bstate%swbotb))
    end subroutine boundary_state_to_variables

    ! ==========================================================================
    ! Solute State Synchronization
    ! ==========================================================================
    subroutine solute_state_from_variables(solu, n_nod, n_lay, n_lev)
        type(solute_state_t), intent(inout) :: solu
        integer, intent(in) :: n_nod, n_lay, n_lev
        integer :: i
        
        ! Configuration switches
        solu%swsolu = swsolu
        solu%swsp = swsp
        solu%swbr = swbr
        solu%swbotbc = swbotbc
        solu%nconc = nconc
        
        ! Concentrations - arrays
        if (allocated(solu%cml) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(solu%cml))
                solu%cml(i) = cml(i)
                solu%cmsy(i) = cmsy(i)
            end do
        end if
        
        ! Concentrations - scalars
        solu%cpond = cpond
        solu%csurf = csurf
        solu%cdrain = cdrain
        solu%cseep = cseep
        solu%cpre = cpre
        solu%cirr = cirr
        solu%cref = cref
        
        ! Cumulative amounts
        solu%sampro = sampro
        solu%samini = samini
        solu%sqbot = sqbot
        solu%sqdra = sqdra
        solu%sqprec = sqprec
        solu%sqirrig = sqirrig
        solu%sqsur = sqsur
        solu%sqrap = sqrap
        solu%dectot = dectot
        solu%rottot = rottot
        solu%solbal = solbal
        
        ! Intermediate amounts
        solu%imsqbot = imsqbot
        solu%imsqdra = imsqdra
        solu%imsqprec = imsqprec
        solu%imsqirrig = imsqirrig
        solu%imdectot = imdectot
        solu%imrottot = imrottot
        solu%isqbot = isqbot
        solu%isqtop = isqtop
        
        ! Macropore solute
        solu%samcra = samcra
        
        ! Age tracer
        solu%AgeGwl1m = AgeGwl1m
        solu%icAgeBot = icAgeBot
        solu%icAgeRot = icAgeRot
        solu%icAgeSur = icAgeSur
        if (allocated(solu%icAgeDra) .and. n_lev > 0) then
            do i = 1, min(n_lev, size(solu%icAgeDra))
                solu%icAgeDra(i) = icAgeDra(i)
            end do
        end if
        
        ! Transport parameters
        solu%ddif = ddif
        solu%frexp = frexp
        solu%tscf = tscf
        solu%dtsolu = dtsolu
        
        ! Decomposition parameters
        solu%gampar = gampar
        solu%bexp = bexp
        solu%rtheta = rtheta
        solu%decsat = decsat
        
        ! Aquifer parameters
        solu%daquif = daquif
        solu%poros = poros
        solu%kfsat = kfsat
        
        ! Salt stress parameters
        solu%salthead = salthead
        solu%saltmax = saltmax
        solu%saltslope = saltslope
        
        ! Parameters per layer
        if (allocated(solu%ldis) .and. n_lay > 0) then
            do i = 1, min(n_lay, size(solu%ldis))
                solu%ldis(i) = ldis(i)
                solu%kf(i) = kf(i)
                solu%decpot(i) = decpot(i)
                solu%fdepth(i) = fdepth(i)
            end do
        end if
        
        ! Tables
        if (allocated(solu%cseeptab)) then
            do i = 1, min(size(solu%cseeptab), size(cseeptab))
                solu%cseeptab(i) = cseeptab(i)
            end do
        end if
        if (allocated(solu%zc) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(solu%zc))
                solu%zc(i) = zc(i)
            end do
        end if
        
        ! Flags
        solu%flsolute = flsolute
        solu%flAgeTracer = flAgeTracer
        
        ! Age tracer boundary/pond state
        solu%Ageirr = Ageirr
        solu%Agedrain = Agedrain
        solu%Agepre = Agepre
        solu%Agepond = Agepond
        solu%Agepondm1 = Agepondm1
        solu%icAgetopupw = icAgetopupw
        solu%icAgetopdwn = icAgetopdwn
        solu%ArMpSs = ArMpSs
        
        call log_debug('sync', 'solute_state_from_variables: swsolu=' // to_str(solu%swsolu))
    end subroutine solute_state_from_variables
    
    subroutine solute_state_to_variables(solu, n_nod, n_lay, n_lev)
        type(solute_state_t), intent(in) :: solu
        integer, intent(in) :: n_nod, n_lay, n_lev
        integer :: i
        
        ! Configuration switches
        swsolu = solu%swsolu
        swsp = solu%swsp
        swbr = solu%swbr
        swbotbc = solu%swbotbc
        nconc = solu%nconc
        
        ! Concentrations - arrays
        if (allocated(solu%cml) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(solu%cml))
                cml(i) = solu%cml(i)
                cmsy(i) = solu%cmsy(i)
            end do
        end if
        
        ! Concentrations - scalars
        cpond = solu%cpond
        csurf = solu%csurf
        cdrain = solu%cdrain
        cseep = solu%cseep
        cpre = solu%cpre
        cirr = solu%cirr
        cref = solu%cref
        
        ! Cumulative amounts
        sampro = solu%sampro
        samini = solu%samini
        sqbot = solu%sqbot
        sqdra = solu%sqdra
        sqprec = solu%sqprec
        sqirrig = solu%sqirrig
        sqsur = solu%sqsur
        sqrap = solu%sqrap
        dectot = solu%dectot
        rottot = solu%rottot
        solbal = solu%solbal
        
        ! Intermediate amounts
        imsqbot = solu%imsqbot
        imsqdra = solu%imsqdra
        imsqprec = solu%imsqprec
        imsqirrig = solu%imsqirrig
        imdectot = solu%imdectot
        imrottot = solu%imrottot
        isqbot = solu%isqbot
        isqtop = solu%isqtop
        
        ! Macropore solute
        samcra = solu%samcra
        
        ! Age tracer
        AgeGwl1m = solu%AgeGwl1m
        icAgeBot = solu%icAgeBot
        icAgeRot = solu%icAgeRot
        icAgeSur = solu%icAgeSur
        if (allocated(solu%icAgeDra) .and. n_lev > 0) then
            do i = 1, min(n_lev, size(solu%icAgeDra))
                icAgeDra(i) = solu%icAgeDra(i)
            end do
        end if
        
        ! Transport parameters
        ddif = solu%ddif
        frexp = solu%frexp
        tscf = solu%tscf
        dtsolu = solu%dtsolu
        
        ! Decomposition parameters
        gampar = solu%gampar
        bexp = solu%bexp
        rtheta = solu%rtheta
        decsat = solu%decsat
        
        ! Aquifer parameters
        daquif = solu%daquif
        poros = solu%poros
        kfsat = solu%kfsat
        
        ! Salt stress parameters
        salthead = solu%salthead
        saltmax = solu%saltmax
        saltslope = solu%saltslope
        
        ! Parameters per layer
        if (allocated(solu%ldis) .and. n_lay > 0) then
            do i = 1, min(n_lay, size(solu%ldis))
                ldis(i) = solu%ldis(i)
                kf(i) = solu%kf(i)
                decpot(i) = solu%decpot(i)
                fdepth(i) = solu%fdepth(i)
            end do
        end if
        
        ! Tables
        if (allocated(solu%cseeptab)) then
            do i = 1, min(size(solu%cseeptab), size(cseeptab))
                cseeptab(i) = solu%cseeptab(i)
            end do
        end if
        if (allocated(solu%zc) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(solu%zc))
                zc(i) = solu%zc(i)
            end do
        end if
        
        ! Flags
        flsolute = solu%flsolute
        flAgeTracer = solu%flAgeTracer
        
        ! Age tracer boundary/pond state
        Ageirr = solu%Ageirr
        Agedrain = solu%Agedrain
        Agepre = solu%Agepre
        Agepond = solu%Agepond
        Agepondm1 = solu%Agepondm1
        icAgetopupw = solu%icAgetopupw
        icAgetopdwn = solu%icAgetopdwn
        ArMpSs = solu%ArMpSs
        
        call log_debug('sync', 'solute_state_to_variables: swsolu=' // to_str(solu%swsolu))
    end subroutine solute_state_to_variables

    ! ==========================================================================
    ! Heat State Synchronization
    ! ==========================================================================
    subroutine heat_state_from_variables(hstate, n_nod, n_lay)
        type(heat_state_t), intent(inout) :: hstate
        integer, intent(in) :: n_nod, n_lay
        integer :: i
        
        ! Configuration switches
        hstate%swhea = swhea
        hstate%swcalt = swcalt
        hstate%swtopbhea = swtopbhea
        hstate%swbotbhea = swbotbhea
        hstate%swfrost = swfrost
        hstate%nheat = nheat
        
        ! Soil temperatures
        if (allocated(hstate%tsoil) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%tsoil))
                hstate%tsoil(i) = tsoil(i)
            end do
        end if
        hstate%tetop = tetop
        hstate%tebot = tebot
        
        ! Thermal properties per compartment
        if (allocated(hstate%heacap) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%heacap))
                hstate%heacap(i) = heacap(i)
                hstate%heacon(i) = heacon(i)
            end do
        end if
        
        ! Frost reduction factor
        if (allocated(hstate%rfcp) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%rfcp))
                hstate%rfcp(i) = rfcp(i)
            end do
        end if
        
        ! Soil composition per compartment
        if (allocated(hstate%fclay) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%fclay))
                hstate%fclay(i) = fclay(i)
                hstate%forg(i) = forg(i)
                hstate%fquartz(i) = fquartz(i)
            end do
        end if
        
        ! Soil composition per layer
        if (allocated(hstate%pclay) .and. n_lay > 0) then
            do i = 1, min(n_lay, size(hstate%pclay))
                hstate%pclay(i) = pclay(i)
                hstate%psand(i) = psand(i)
                hstate%psilt(i) = psilt(i)
                hstate%orgmat(i) = orgmat(i)
            end do
        end if
        
        ! Boundary conditions
        hstate%tmean = tmean
        hstate%tampli = tampli
        hstate%timref = timref
        hstate%ddamp = ddamp
        
        ! Boundary condition tables
        if (allocated(hstate%tembtab)) then
            do i = 1, min(size(hstate%tembtab), size(tembtab))
                hstate%tembtab(i) = tembtab(i)
            end do
        end if
        if (allocated(hstate%temtoptab)) then
            do i = 1, min(size(hstate%temtoptab), size(temtoptab))
                hstate%temtoptab(i) = temtoptab(i)
            end do
        end if
        if (allocated(hstate%zh) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%zh))
                hstate%zh(i) = zh(i)
            end do
        end if
        
        ! Frost
        hstate%zfrosttop = zfrosttop
        hstate%zfrostbot = zfrostbot
        hstate%tfroststa = tfroststa
        hstate%tfrostend = tfrostend
        hstate%nodfrostbot = nodfrostbot
        
        ! Flags
        hstate%fltemperature = fltemperature
        
        call log_debug('sync', 'heat_state_from_variables: swhea=' // to_str(hstate%swhea))
    end subroutine heat_state_from_variables
    
    subroutine heat_state_to_variables(hstate, n_nod, n_lay)
        type(heat_state_t), intent(in) :: hstate
        integer, intent(in) :: n_nod, n_lay
        integer :: i
        
        ! Configuration switches
        swhea = hstate%swhea
        swcalt = hstate%swcalt
        swtopbhea = hstate%swtopbhea
        swbotbhea = hstate%swbotbhea
        swfrost = hstate%swfrost
        nheat = hstate%nheat
        
        ! Soil temperatures
        if (allocated(hstate%tsoil) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%tsoil))
                tsoil(i) = hstate%tsoil(i)
            end do
        end if
        tetop = hstate%tetop
        tebot = hstate%tebot
        
        ! Thermal properties per compartment
        if (allocated(hstate%heacap) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%heacap))
                heacap(i) = hstate%heacap(i)
                heacon(i) = hstate%heacon(i)
            end do
        end if
        
        ! Frost reduction factor
        if (allocated(hstate%rfcp) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%rfcp))
                rfcp(i) = hstate%rfcp(i)
            end do
        end if
        
        ! Soil composition per compartment
        if (allocated(hstate%fclay) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%fclay))
                fclay(i) = hstate%fclay(i)
                forg(i) = hstate%forg(i)
                fquartz(i) = hstate%fquartz(i)
            end do
        end if
        
        ! Soil composition per layer
        if (allocated(hstate%pclay) .and. n_lay > 0) then
            do i = 1, min(n_lay, size(hstate%pclay))
                pclay(i) = hstate%pclay(i)
                psand(i) = hstate%psand(i)
                psilt(i) = hstate%psilt(i)
                orgmat(i) = hstate%orgmat(i)
            end do
        end if
        
        ! Boundary conditions
        tmean = hstate%tmean
        tampli = hstate%tampli
        timref = hstate%timref
        ddamp = hstate%ddamp
        
        ! Boundary condition tables
        if (allocated(hstate%tembtab)) then
            do i = 1, min(size(hstate%tembtab), size(tembtab))
                tembtab(i) = hstate%tembtab(i)
            end do
        end if
        if (allocated(hstate%temtoptab)) then
            do i = 1, min(size(hstate%temtoptab), size(temtoptab))
                temtoptab(i) = hstate%temtoptab(i)
            end do
        end if
        if (allocated(hstate%zh) .and. n_nod > 0) then
            do i = 1, min(n_nod, size(hstate%zh))
                zh(i) = hstate%zh(i)
            end do
        end if
        
        ! Frost
        zfrosttop = hstate%zfrosttop
        zfrostbot = hstate%zfrostbot
        tfroststa = hstate%tfroststa
        tfrostend = hstate%tfrostend
        nodfrostbot = hstate%nodfrostbot
        
        ! Flags
        fltemperature = hstate%fltemperature
        
        call log_debug('sync', 'heat_state_to_variables: swhea=' // to_str(hstate%swhea))
    end subroutine heat_state_to_variables

end module swap_state_sync
