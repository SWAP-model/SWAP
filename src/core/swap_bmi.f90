! swap_bmi_refactored.f90
! Basic Model Interface (BMI) implementation for SWAP with state object
! 
! This is a PROOF-OF-CONCEPT refactoring showing state object usage
! while maintaining backward compatibility

module swap_bmi
    use iso_c_binding, only: c_int, c_char, c_double, c_ptr, c_loc, c_null_char
    use variables
    use swap_exchange
    use swap_state  ! NEW: State object module
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
    
    ! NEW: State object (for now, single global instance for backward compat)
    type(swap_state_t), save :: global_state
    
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
        
        ! NEW: Initialize state object
        global_state = swap_state_create()
        
        ! Convert C string to Fortran string
        length = 0
        do i = 1, 256
            if (config_file(i) == c_null_char) exit
            fortran_config_file(i:i) = config_file(i)
            length = i
        end do
        
        ! Call SWAP initialization (iTask = 1)
        call swap(global_state%iCaller, 1)
        
        ! NEW: Sync state from global variables
        call swap_state_sync_from_globals(global_state)
        
        global_state%initialized = .true.
        status = BMI_SUCCESS
    end function bmi_initialize
    
    !-------------------------------------------------------------------
    ! Initialize without reading files (for batch processing)
    !-------------------------------------------------------------------
    function bmi_initialize_memory() result(status) bind(C, name="initialize_memory")
        integer(c_int) :: status
        
        ! NEW: Initialize state object
        global_state = swap_state_create()
        
        ! Initialize minimal SWAP state without file reading
        call initialize()
        
        global_state%memory_mode = .true.
        global_state%initialized = .true.
        status = BMI_SUCCESS
    end function bmi_initialize_memory
    
    !-------------------------------------------------------------------
    ! Run one time step
    !-------------------------------------------------------------------
    function bmi_update() result(status) bind(C, name="update")
        integer(c_int) :: status
        
        if (.not. global_state%initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Run SWAP dynamic loop (iTask = 2)
        call swap(global_state%iCaller, 2)
        
        ! NEW: Sync state from global variables
        call swap_state_sync_from_globals(global_state)
        
        status = BMI_SUCCESS
    end function bmi_update
    
    !-------------------------------------------------------------------
    ! Run one day (more practical than update for daily models)
    !-------------------------------------------------------------------
    function bmi_update_day() result(status) bind(C, name="update_day")
        integer(c_int) :: status
        
        if (.not. global_state%initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Execute one day of simulation
        flrunend = .false.
        fldayend = .false.
        
        ! Run until end of day
        do while (.not. fldayend)
            call swap(global_state%iCaller, 2)
        end do
        
        ! NEW: Sync state from global variables
        call swap_state_sync_from_globals(global_state)
        
        status = BMI_SUCCESS
    end function bmi_update_day
    
    !-------------------------------------------------------------------
    ! Finalize the model
    !-------------------------------------------------------------------
    function bmi_finalize() result(status) bind(C, name="finalize")
        integer(c_int) :: status
        
        if (.not. global_state%initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ! Call SWAP closure (iTask = 3)
        call swap(global_state%iCaller, 3)
        
        ! NEW: Clean up state object
        call swap_state_destroy(global_state)
        global_state%initialized = .false.
        
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
        
        ! NEW: Return from state object (after sync)
        call swap_state_sync_from_globals(global_state)
        time = global_state%t1900
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
        
        character(len=*), parameter :: component_name = "SWAP 4.2.0 (State Object POC)"
        integer :: i
        
        do i = 1, len_trim(component_name)
            name(i) = component_name(i:i)
        end do
        name(len_trim(component_name)+1) = c_null_char
        
        status = BMI_SUCCESS
    end function bmi_get_component_name

!=======================================================================
! BMI GETTER FUNCTIONS
!=======================================================================

    function bmi_get_value_ptr_theta(ptr) result(status) bind(C, name="get_value_ptr_theta")
        type(c_ptr), intent(out) :: ptr
        integer(c_int) :: status
        
        if (.not. global_state%initialized) then
            status = BMI_FAILURE
            return
        end if
        
        ptr = c_loc(theta(1))
        status = BMI_SUCCESS
    end function bmi_get_value_ptr_theta
    
    function bmi_get_value_gwl(gwl_out) result(status) bind(C, name="get_value_gwl")
        real(c_double), intent(out) :: gwl_out
        integer(c_int) :: status
        
        gwl_out = gwl
        status = BMI_SUCCESS
    end function bmi_get_value_gwl
    
    function bmi_get_grid_size(grid_size) result(status) bind(C, name="get_grid_size")
        integer(c_int), intent(out) :: grid_size
        integer(c_int) :: status
        
        grid_size = numnod
        status = BMI_SUCCESS
    end function bmi_get_grid_size
    
    !-------------------------------------------------------------------
    ! Get scalar: cumulative gross rainfall (cm)
    !-------------------------------------------------------------------
    function bmi_get_value_cgrai(cgrai_out) result(status) bind(C, name="get_value_cgrai")
        real(c_double), intent(out) :: cgrai_out
        integer(c_int) :: status
        
        cgrai_out = cgrai
        status = BMI_SUCCESS
    end function bmi_get_value_cgrai
    
    !-------------------------------------------------------------------
    ! Get scalar: cumulative net rainfall (cm)
    !-------------------------------------------------------------------
    function bmi_get_value_cnrai(cnrai_out) result(status) bind(C, name="get_value_cnrai")
        real(c_double), intent(out) :: cnrai_out
        integer(c_int) :: status
        
        cnrai_out = cnrai
        status = BMI_SUCCESS
    end function bmi_get_value_cnrai
    
    !-------------------------------------------------------------------
    ! Get scalar: current year
    !-------------------------------------------------------------------
    function bmi_get_value_year(year_out) result(status) bind(C, name="get_value_year")
        integer(c_int), intent(out) :: year_out
        integer(c_int) :: status
        
        year_out = iyear
        status = BMI_SUCCESS
    end function bmi_get_value_year
    
    !-------------------------------------------------------------------
    ! Get scalar: current day of year
    !-------------------------------------------------------------------
    function bmi_get_value_daynr(daynr_out) result(status) bind(C, name="get_value_daynr")
        integer(c_int), intent(out) :: daynr_out
        integer(c_int) :: status
        
        daynr_out = daynr
        status = BMI_SUCCESS
    end function bmi_get_value_daynr

end module swap_bmi
