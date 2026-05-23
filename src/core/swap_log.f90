! ==============================================================================
! SWAP Logging Module
! ==============================================================================
! Simple logging facility for debugging SWAP model execution.
!
! Features:
!   - Multiple log levels (DEBUG, INFO, WARN, ERROR)
!   - Optional file output
!   - Timestamps
!   - Module/subroutine context
!   - Compile-time disabling via preprocessor
!
! Usage:
!   use swap_log
!   call log_init(log_level=LOG_DEBUG, log_file='swap_debug.log')
!   call log_debug('state', 'Initializing soil state with n=' // to_str(n))
!   call log_info('state', 'State initialized successfully')
!   call log_warn('solver', 'Convergence slow, iter=' // to_str(iter))
!   call log_error('io', 'Failed to open file: ' // trim(filename))
!   call log_close()
!
! Author: SWAP Development Team
! Date: 2026-01-31
! ==============================================================================

module swap_log
    implicit none
    private
    
    ! Log levels (lower number = more verbose)
    integer, parameter, public :: LOGLEVEL_DEBUG = 10
    integer, parameter, public :: LOGLEVEL_INFO  = 20
    integer, parameter, public :: LOGLEVEL_WARN  = 30
    integer, parameter, public :: LOGLEVEL_ERROR = 40
    integer, parameter, public :: LOGLEVEL_NONE  = 100
    
    ! Module state
    integer :: current_level = LOGLEVEL_INFO
    integer :: log_unit = -1
    logical :: log_to_file = .false.
    logical :: log_to_stdout = .true.
    logical :: log_initialized = .false.
    logical :: include_timestamp = .true.
    
    ! Public interface
    public :: log_init
    public :: log_close
    public :: log_debug
    public :: log_info
    public :: log_warn
    public :: log_error
    public :: log_set_level
    public :: log_message
    public :: log_unit_handle
    public :: to_str
    
    ! Generic interface for to_str
    interface to_str
        module procedure int_to_str
        module procedure real_to_str
        module procedure real8_to_str
        module procedure logical_to_str
    end interface to_str

contains

    ! ===========================================================================
    ! Initialization and Cleanup
    ! ===========================================================================
    
    subroutine log_init(log_level, log_file, to_stdout, timestamps)
        !> Initialize the logging system
        integer, intent(in), optional :: log_level        ! Minimum level to log
        character(len=*), intent(in), optional :: log_file ! File to write logs to
        logical, intent(in), optional :: to_stdout        ! Also write to stdout
        logical, intent(in), optional :: timestamps       ! Include timestamps
        
        integer :: ios
        
        ! Set log level
        if (present(log_level)) then
            current_level = log_level
        else
            current_level = LOGLEVEL_INFO
        end if
        
        ! Set stdout option
        if (present(to_stdout)) then
            log_to_stdout = to_stdout
        else
            log_to_stdout = .true.
        end if
        
        ! Set timestamp option
        if (present(timestamps)) then
            include_timestamp = timestamps
        else
            include_timestamp = .true.
        end if
        
        ! Open log file if specified
        if (present(log_file)) then
            open(newunit=log_unit, file=log_file, status='replace', &
                 action='write', iostat=ios)
            if (ios == 0) then
                log_to_file = .true.
            else
                log_to_file = .false.
                log_unit = -1
                write(*,'(A)') 'WARNING: Could not open log file: ' // trim(log_file)
            end if
        end if
        
        log_initialized = .true.
        
        ! Log initialization message
        call log_info('log', 'Logging initialized at level ' // level_name(current_level))
        
    end subroutine log_init
    
    subroutine log_close()
        !> Close the logging system and flush buffers
        if (log_to_file .and. log_unit /= -1) then
            call log_info('log', 'Logging closed')
            close(log_unit)
            log_unit = -1
            log_to_file = .false.
        end if
        log_initialized = .false.
    end subroutine log_close
    
    function log_unit_handle() result(unit)
        !> Return the raw file unit number swap_log uses for the log file.
        !! Use only when interop with legacy code that calls write(unit, ...)
        !! directly is unavoidable (e.g. msw1eic's quarantined error path).
        !! Prefer log_debug/info/warn/error for new code.
        integer :: unit
        unit = log_unit
    end function log_unit_handle

    subroutine log_set_level(level)
        !> Change the current log level
        integer, intent(in) :: level
        current_level = level
    end subroutine log_set_level

    ! ===========================================================================
    ! Logging Functions
    ! ===========================================================================
    
    subroutine log_debug(context, message)
        !> Log a debug message
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: message
        call log_message(LOGLEVEL_DEBUG, context, message)
    end subroutine log_debug
    
    subroutine log_info(context, message)
        !> Log an info message
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: message
        call log_message(LOGLEVEL_INFO, context, message)
    end subroutine log_info
    
    subroutine log_warn(context, message)
        !> Log a warning message
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: message
        call log_message(LOGLEVEL_WARN, context, message)
    end subroutine log_warn
    
    subroutine log_error(context, message)
        !> Log an error message
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: message
        call log_message(LOGLEVEL_ERROR, context, message)
    end subroutine log_error
    
    subroutine log_message(level, context, message)
        !> Core logging routine
        integer, intent(in) :: level
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: message
        
        character(len=512) :: formatted_msg
        character(len=23) :: timestamp_str
        character(len=5) :: level_str
        
        ! Check if we should log this message
        if (level < current_level) return
        
        ! Format level string
        level_str = level_name(level)
        
        ! Build formatted message
        if (include_timestamp) then
            timestamp_str = get_timestamp()
            write(formatted_msg, '(A,1X,A5,1X,A,": ",A)') &
                trim(timestamp_str), level_str, trim(context), trim(message)
        else
            write(formatted_msg, '(A5,1X,A,": ",A)') &
                level_str, trim(context), trim(message)
        end if
        
        ! Output to stdout
        if (log_to_stdout) then
            write(*,'(A)') trim(formatted_msg)
        end if
        
        ! Output to file
        if (log_to_file .and. log_unit /= -1) then
            write(log_unit,'(A)') trim(formatted_msg)
            flush(log_unit)  ! Ensure immediate write for debugging
        end if
        
    end subroutine log_message

    ! ===========================================================================
    ! Helper Functions
    ! ===========================================================================
    
    function level_name(level) result(name)
        !> Get human-readable name for log level
        integer, intent(in) :: level
        character(len=5) :: name
        
        select case (level)
            case (LOGLEVEL_DEBUG)
                name = 'DEBUG'
            case (LOGLEVEL_INFO)
                name = 'INFO '
            case (LOGLEVEL_WARN)
                name = 'WARN '
            case (LOGLEVEL_ERROR)
                name = 'ERROR'
            case default
                name = '?????'
        end select
    end function level_name
    
    function get_timestamp() result(ts)
        !> Get current timestamp as string
        character(len=23) :: ts
        character(len=8) :: date_str
        character(len=10) :: time_str
        
        call date_and_time(date=date_str, time=time_str)
        
        ! Format: YYYY-MM-DD HH:MM:SS.sss
        write(ts, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2,".",A3)') &
            date_str(1:4), date_str(5:6), date_str(7:8), &
            time_str(1:2), time_str(3:4), time_str(5:6), time_str(8:10)
    end function get_timestamp
    
    function int_to_str(val) result(str)
        !> Convert integer to string
        integer, intent(in) :: val
        character(len=20) :: str
        write(str, '(I0)') val
        str = adjustl(str)
    end function int_to_str
    
    function real_to_str(val) result(str)
        !> Convert real to string
        real, intent(in) :: val
        character(len=20) :: str
        write(str, '(G12.5)') val
        str = adjustl(str)
    end function real_to_str
    
    function real8_to_str(val) result(str)
        !> Convert real(8) to string
        real(8), intent(in) :: val
        character(len=24) :: str
        write(str, '(G15.7)') val
        str = adjustl(str)
    end function real8_to_str
    
    function logical_to_str(val) result(str)
        !> Convert logical to string
        logical, intent(in) :: val
        character(len=5) :: str
        if (val) then
            str = 'true'
        else
            str = 'false'
        end if
    end function logical_to_str

end module swap_log
