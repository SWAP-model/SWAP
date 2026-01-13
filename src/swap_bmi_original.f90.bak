! swap_bmi.f90
! Basic Model Interface (BMI) implementation for SWAP
! 
! This module provides a standard interface for running SWAP as a library
! from Python or other programming languages. It follows the CSDMS BMI
! specification: https://bmi.csdms.io/
!
! Author: Zawadzki M.
! Date: 2026-01-08
!
! Primary use cases:
!   1. Run SWAP from Python without file I/O overhead
!   2. Parameter estimation and calibration
!   3. Ensemble simulations
!   4. Quasi-spatial analysis
!   5. Model coupling (e.g., with surface water models)

module swap_bmi
    use iso_c_binding, only: c_int, c_char, c_double, c_ptr, c_loc, c_null_char
    use variables
    use swap_exchange
    implicit none
    
    ! BMI status codes
    integer(c_int), parameter :: BMI_SUCCESS = 0
    integer(c_int), parameter :: BMI_FAILURE = 1
    
    ! BMI string length parameters (CSDMS standard)
    integer(c_int), parameter :: BMI_MAX_COMPONENT_NAME = 2048
    integer(c_int), parameter :: BMI_MAX_VAR_NAME = 2048
    integer(c_int), parameter :: BMI_MAX_TYPE_NAME = 2048
    integer(c_int), parameter :: BMI_MAX_UNITS_NAME = 2048
    
    ! Expose these for Python bindings
    integer(c_int), bind(C, name="BMI_LENCOMPONENTNAME") :: BMI_LENCOMPONENTNAME = 2048
    integer(c_int), bind(C, name="BMI_LENVARNAME") :: BMI_LENVARNAME = 2048
    integer(c_int), bind(C, name="BMI_LENVARTYPE") :: BMI_LENVARTYPE = 2048
    integer(c_int), bind(C, name="BMI_LENVARUNITS") :: BMI_LENVARUNITS = 2048
    
    ! Internal state
    logical, save :: bmi_initialized = .false.
    logical, save :: bmi_memory_mode = .false.  ! Skip file I/O when true
    integer, save :: bmi_iCaller = 0
    
contains

!=======================================================================
! BMI LIFECYCLE FUNCTIONS
!=======================================================================

    !-------------------------------------------------------------------
    ! Initialize the model
    !-------------------------------------------------------------------
    function bmi_initialize(config_file) result(status) bind(C, name="initialize")
        character(kind=c_char), intent(in) :: config_file(*)
        integer(c_int) :: status
        
        character(len=256) :: fortran_config_file
        integer :: i, length
        
        ! Convert C string to Fortran string
        length = 0
        do i = 1, 256
            if (config_file(i) == c_null_char) exit
            fortran_config_file(i:i) = config_file(i)
            length = i
        end do
        
        ! Call SWAP initialization (iTask = 1)
        call swap(bmi_iCaller, 1)
        
        bmi_initialized = .true.
        status = BMI_SUCCESS
    end function bmi_initialize
    
    !-------------------------------------------------------------------
    ! Initialize without reading files (for batch processing)
    !-------------------------------------------------------------------
    function bmi_initialize_memory() result(status) bind(C, name="initialize_memory")
        integer(c_int) :: status
        
        ! Initialize minimal SWAP state without file reading
        call initialize()
        
        bmi_memory_mode = .true.
        bmi_initialized = .true.
        status = BMI_SUCCESS
    end function bmi_initialize_memory
    
    !-------------------------------------------------------------------
    ! Run one time step
    !-------------------------------------------------------------------
    function bmi_update() result(status) bind(C, name="update")
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Run SWAP dynamic loop (iTask = 2)
        ! This will run until flrunend becomes true
        call swap(bmi_iCaller, 2)
        
        status = BMI_SUCCESS
    end function bmi_update
    
    !-------------------------------------------------------------------
    ! Run one day (more practical than update for daily models)
    !-------------------------------------------------------------------
    function bmi_update_day() result(status) bind(C, name="update_day")
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Execute one day of simulation
        ! Set flrunend to stop after one day
        flrunend = .false.
        fldayend = .false.
        
        ! Run until end of day
        do while (.not. fldayend)
            call swap(bmi_iCaller, 2)
        end do
        
        status = BMI_SUCCESS
    end function bmi_update_day
    
    !-------------------------------------------------------------------
    ! Finalize the model
    !-------------------------------------------------------------------
    function bmi_finalize() result(status) bind(C, name="finalize")
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Call SWAP closure (iTask = 3)
        call swap(bmi_iCaller, 3)
        
        bmi_initialized = .false.
        status = BMI_SUCCESS
    end function bmi_finalize

!=======================================================================
! BMI TIME FUNCTIONS
!=======================================================================

    !-------------------------------------------------------------------
    ! Get current model time (days since 1900)
    !-------------------------------------------------------------------
    function bmi_get_current_time(time) result(status) bind(C, name="get_current_time")
        real(c_double), intent(out) :: time
        integer(c_int) :: status
        
        time = t1900
        status = BMI_SUCCESS
    end function bmi_get_current_time
    
    !-------------------------------------------------------------------
    ! Get model start time
    !-------------------------------------------------------------------
    function bmi_get_start_time(time) result(status) bind(C, name="get_start_time")
        real(c_double), intent(out) :: time
        integer(c_int) :: status
        
        time = tstart
        status = BMI_SUCCESS
    end function bmi_get_start_time
    
    !-------------------------------------------------------------------
    ! Get model end time
    !-------------------------------------------------------------------
    function bmi_get_end_time(time) result(status) bind(C, name="get_end_time")
        real(c_double), intent(out) :: time
        integer(c_int) :: status
        
        time = tend
        status = BMI_SUCCESS
    end function bmi_get_end_time
    
    !-------------------------------------------------------------------
    ! Get current time step
    !-------------------------------------------------------------------
    function bmi_get_time_step(time_step) result(status) bind(C, name="get_time_step")
        real(c_double), intent(out) :: time_step
        integer(c_int) :: status
        
        time_step = dt
        status = BMI_SUCCESS
    end function bmi_get_time_step
    
    !-------------------------------------------------------------------
    ! Get time units
    !-------------------------------------------------------------------
    function bmi_get_time_units(units) result(status) bind(C, name="get_time_units")
        character(kind=c_char), intent(out) :: units(BMI_MAX_UNITS_NAME)
        integer(c_int) :: status
        
        character(len=*), parameter :: units_str = "days since 1900-01-01"
        integer :: i
        
        do i = 1, len_trim(units_str)
            units(i) = units_str(i:i)
        end do
        units(len_trim(units_str)+1) = c_null_char
        
        status = BMI_SUCCESS
    end function bmi_get_time_units

!=======================================================================
! BMI VARIABLE INFORMATION
!=======================================================================

    !-------------------------------------------------------------------
    ! Get component name
    !-------------------------------------------------------------------
    function bmi_get_component_name(name) result(status) bind(C, name="get_component_name")
        character(kind=c_char), intent(out) :: name(BMI_MAX_COMPONENT_NAME)
        integer(c_int) :: status
        
        character(len=*), parameter :: component_name = "SWAP 4.2.0"
        integer :: i
        
        do i = 1, len_trim(component_name)
            name(i) = component_name(i:i)
        end do
        name(len_trim(component_name)+1) = c_null_char
        
        status = BMI_SUCCESS
    end function bmi_get_component_name

!=======================================================================
! BMI GETTER FUNCTIONS - Direct memory access for efficiency
!=======================================================================

    !-------------------------------------------------------------------
    ! Get pointer to soil moisture (theta) array
    !-------------------------------------------------------------------
    function bmi_get_value_ptr_theta(ptr) result(status) bind(C, name="get_value_ptr_theta")
        type(c_ptr), intent(out) :: ptr
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Return pointer to theta array (volumetric water content)
        ptr = c_loc(theta(1))
        status = BMI_SUCCESS
    end function bmi_get_value_ptr_theta
    
    !-------------------------------------------------------------------
    ! Get pointer to pressure head (h) array
    !-------------------------------------------------------------------
    function bmi_get_value_ptr_h(ptr) result(status) bind(C, name="get_value_ptr_h")
        type(c_ptr), intent(out) :: ptr
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ptr = c_loc(h(1))
        status = BMI_SUCCESS
    end function bmi_get_value_ptr_h
    
    !-------------------------------------------------------------------
    ! Get pointer to layer thickness (dz) array
    !-------------------------------------------------------------------
    function bmi_get_value_ptr_dz(ptr) result(status) bind(C, name="get_value_ptr_dz")
        type(c_ptr), intent(out) :: ptr
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ptr = c_loc(dz(1))
        status = BMI_SUCCESS
    end function bmi_get_value_ptr_dz
    
    !-------------------------------------------------------------------
    ! Get pointer to depth (z) array
    !-------------------------------------------------------------------
    function bmi_get_value_ptr_z(ptr) result(status) bind(C, name="get_value_ptr_z")
        type(c_ptr), intent(out) :: ptr
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ptr = c_loc(z(1))
        status = BMI_SUCCESS
    end function bmi_get_value_ptr_z
    
    !-------------------------------------------------------------------
    ! Get pointer to root water uptake array
    !-------------------------------------------------------------------
    function bmi_get_value_ptr_rwu(ptr) result(status) bind(C, name="get_value_ptr_rwu")
        type(c_ptr), intent(out) :: ptr
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ptr = c_loc(inqrot(1))
        status = BMI_SUCCESS
    end function bmi_get_value_ptr_rwu
    
    !-------------------------------------------------------------------
    ! Get scalar: groundwater level
    !-------------------------------------------------------------------
    function bmi_get_value_gwl(gwl_out) result(status) bind(C, name="get_value_gwl")
        real(c_double), intent(out) :: gwl_out
        integer(c_int) :: status
        
        gwl_out = gwl
        status = BMI_SUCCESS
    end function bmi_get_value_gwl
    
    !-------------------------------------------------------------------
    ! Get scalar: number of nodes
    !-------------------------------------------------------------------
    function bmi_get_grid_size(grid_size) result(status) bind(C, name="get_grid_size")
        integer(c_int), intent(out) :: grid_size
        integer(c_int) :: status
        
        grid_size = numnod
        status = BMI_SUCCESS
    end function bmi_get_grid_size
    
    !-------------------------------------------------------------------
    ! Get scalar: actual transpiration
    !-------------------------------------------------------------------
    function bmi_get_value_tact(tact) result(status) bind(C, name="get_value_tact")
        real(c_double), intent(out) :: tact
        integer(c_int) :: status
        
        tact = sum(inqrot(1:numnod))
        status = BMI_SUCCESS
    end function bmi_get_value_tact
    
    !-------------------------------------------------------------------
    ! Get scalar: potential transpiration
    !-------------------------------------------------------------------
    function bmi_get_value_tpot(tpot) result(status) bind(C, name="get_value_tpot")
        real(c_double), intent(out) :: tpot
        integer(c_int) :: status
        
        tpot = iptra
        status = BMI_SUCCESS
    end function bmi_get_value_tpot

!=======================================================================
! BMI SETTER FUNCTIONS - For parameter setting
!=======================================================================

    !-------------------------------------------------------------------
    ! Set soil hydraulic parameters: Ksat (cm/day)
    !-------------------------------------------------------------------
    function bmi_set_value_ksat(src) result(status) bind(C, name="set_value_ksat")
        real(c_double), intent(in) :: src(*)
        integer(c_int) :: status
        integer :: node
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Set saturated hydraulic conductivity for each node
        ! cofgen(3, node) is Ksat in SWAP's internal structure
        do node = 1, numnod
            cofgen(3, node) = src(node)
        end do
        
        status = BMI_SUCCESS
    end function bmi_set_value_ksat
    
    !-------------------------------------------------------------------
    ! Set Van Genuchten alpha parameter (1/cm)
    !-------------------------------------------------------------------
    function bmi_set_value_alpha(src) result(status) bind(C, name="set_value_alpha")
        real(c_double), intent(in) :: src(*)
        integer(c_int) :: status
        integer :: node
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! cofgen(4, node) is alpha
        do node = 1, numnod
            cofgen(4, node) = src(node)
        end do
        
        status = BMI_SUCCESS
    end function bmi_set_value_alpha
    
    !-------------------------------------------------------------------
    ! Set Van Genuchten n parameter (-)
    !-------------------------------------------------------------------
    function bmi_set_value_n(src) result(status) bind(C, name="set_value_n")
        real(c_double), intent(in) :: src(*)
        integer(c_int) :: status
        integer :: node
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! cofgen(5, node) is n
        do node = 1, numnod
            cofgen(5, node) = src(node)
        end do
        
        status = BMI_SUCCESS
    end function bmi_set_value_n
    
    !-------------------------------------------------------------------
    ! Set residual water content theta_r (-)
    !-------------------------------------------------------------------
    function bmi_set_value_thetar(src) result(status) bind(C, name="set_value_thetar")
        real(c_double), intent(in) :: src(*)
        integer(c_int) :: status
        integer :: node
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! cofgen(1, node) is theta_r
        do node = 1, numnod
            cofgen(1, node) = src(node)
            thetar(node) = src(node)
        end do
        
        status = BMI_SUCCESS
    end function bmi_set_value_thetar
    
    !-------------------------------------------------------------------
    ! Set saturated water content theta_s (-)
    !-------------------------------------------------------------------
    function bmi_set_value_thetas(src) result(status) bind(C, name="set_value_thetas")
        real(c_double), intent(in) :: src(*)
        integer(c_int) :: status
        integer :: node
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! cofgen(2, node) is theta_s
        do node = 1, numnod
            cofgen(2, node) = src(node)
            thetas(node) = src(node)
        end do
        
        status = BMI_SUCCESS
    end function bmi_set_value_thetas
    
    !-------------------------------------------------------------------
    ! Set daily rainfall (mm/day)
    !-------------------------------------------------------------------
    function bmi_set_value_rainfall(rain) result(status) bind(C, name="set_value_rainfall")
        real(c_double), intent(in) :: rain
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Set rainfall for current day
        ! arai contains daily rainfall array
        arai(daynr) = rain
        
        status = BMI_SUCCESS
    end function bmi_set_value_rainfall
    
    !-------------------------------------------------------------------
    ! Set daily ET reference (mm/day)
    !-------------------------------------------------------------------
    function bmi_set_value_etref(etref_val) result(status) bind(C, name="set_value_etref")
        real(c_double), intent(in) :: etref_val
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        aetr(daynr) = etref_val
        
        status = BMI_SUCCESS
    end function bmi_set_value_etref
    
    !-------------------------------------------------------------------
    ! Set daily min temperature (°C)
    !-------------------------------------------------------------------
    function bmi_set_value_tmin(tmin_val) result(status) bind(C, name="set_value_tmin")
        real(c_double), intent(in) :: tmin_val
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        atmn(daynr) = tmin_val
        
        status = BMI_SUCCESS
    end function bmi_set_value_tmin
    
    !-------------------------------------------------------------------
    ! Set daily max temperature (°C)
    !-------------------------------------------------------------------
    function bmi_set_value_tmax(tmax_val) result(status) bind(C, name="set_value_tmax")
        real(c_double), intent(in) :: tmax_val
        integer(c_int) :: status
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        atmx(daynr) = tmax_val
        
        status = BMI_SUCCESS
    end function bmi_set_value_tmax
    
    !-------------------------------------------------------------------
    ! Set initial soil moisture profile
    !-------------------------------------------------------------------
    function bmi_set_value_theta(src) result(status) bind(C, name="set_value_theta")
        real(c_double), intent(in) :: src(*)
        integer(c_int) :: status
        integer :: node
        
        if (.not. bmi_initialized) then
            status = BMI_FAILURE
            return
        end if
        
        do node = 1, numnod
            theta(node) = src(node)
        end do
        
        status = BMI_SUCCESS
    end function bmi_set_value_theta

!=======================================================================
! HELPER FUNCTIONS
!=======================================================================

    !-------------------------------------------------------------------
    ! Get version string
    !-------------------------------------------------------------------
    function bmi_get_version(version) result(status) bind(C, name="get_version")
        character(kind=c_char), intent(out) :: version(BMI_MAX_COMPONENT_NAME)
        integer(c_int) :: status
        
        character(len=*), parameter :: version_str = "SWAP 4.2.0 with BMI"
        integer :: i
        
        do i = 1, len_trim(version_str)
            version(i) = version_str(i:i)
        end do
        version(len_trim(version_str)+1) = c_null_char
        
        status = BMI_SUCCESS
    end function bmi_get_version

end module swap_bmi
