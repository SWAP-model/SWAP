!> Diagnostics configuration layer (Phase A).
!!
!! Pure resolution of logging settings from layered sources with the
!! precedence  C-API > env > TOML > built-in default. This
!! phase owns only the *config* (level/file/routing) + a precedence merge;
!! the per-instance state-borne diagnostics record is a later phase.
module diagnostics_mod
   use swap_log, only: LOGLEVEL_DEBUG, LOGLEVEL_INFO, LOGLEVEL_WARN, &
                       LOGLEVEL_ERROR, LOGLEVEL_NONE
   implicit none
   private

   public :: diagnostics_config_t, diag_overrides_t
   public :: log_level_from_name
   public :: merge_overrides, resolve_diagnostics_config
   public :: read_env_overrides, default_cli_config, default_embedded_config

   !> Fully-resolved logging settings handed to log_init.
   type :: diagnostics_config_t
      integer                       :: level      = LOGLEVEL_INFO
      logical                       :: to_stdout  = .true.
      logical                       :: to_stderr  = .true.
      logical                       :: timestamps = .false.
      character(len=:), allocatable :: log_file       ! unallocated => no file
   end type diagnostics_config_t

   !> A sparse set of overrides from one source. Only fields whose
   !! has_* flag is .true. override the base config in merge_overrides.
   !! The value fields (level/to_stdout/...) are sentinels read only when
   !! their has_* flag is .true.
   type :: diag_overrides_t
      logical                       :: has_level      = .false.
      integer                       :: level          = LOGLEVEL_INFO
      logical                       :: has_stdout     = .false.
      logical                       :: to_stdout      = .true.
      logical                       :: has_stderr     = .false.
      logical                       :: to_stderr      = .true.
      logical                       :: has_timestamps = .false.
      logical                       :: timestamps     = .false.
      logical                       :: has_file       = .false.
      character(len=:), allocatable :: log_file
   end type diag_overrides_t

contains

   !> Map a level name (case-insensitive) to a LOGLEVEL_* value.
   !! Unknown/empty => LOGLEVEL_INFO (safe default).
   pure function log_level_from_name(name) result(level)
      character(len=*), intent(in) :: name
      integer :: level
      select case (upcase(trim(adjustl(name))))
      case ('DEBUG');           level = LOGLEVEL_DEBUG
      case ('INFO');            level = LOGLEVEL_INFO
      case ('WARN', 'WARNING'); level = LOGLEVEL_WARN
      case ('ERROR');           level = LOGLEVEL_ERROR
      case ('NONE', 'OFF');     level = LOGLEVEL_NONE
      case default;             level = LOGLEVEL_INFO
      end select
   end function log_level_from_name

   !> Apply a sparse override set onto a base config (only set fields win).
   pure function merge_overrides(base, ov) result(cfg)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t),     intent(in) :: ov
      type(diagnostics_config_t)             :: cfg
      cfg = base
      if (ov%has_level)      cfg%level      = ov%level
      if (ov%has_stdout)     cfg%to_stdout  = ov%to_stdout
      if (ov%has_stderr)     cfg%to_stderr  = ov%to_stderr
      if (ov%has_timestamps) cfg%timestamps = ov%timestamps
      if (ov%has_file)       cfg%log_file   = ov%log_file
   end function merge_overrides

   !> Resolve the final config with precedence  capi > env > toml > base.
   pure function resolve_diagnostics_config(base, toml, env, capi) result(cfg)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t),     intent(in) :: toml, env, capi
      type(diagnostics_config_t)             :: cfg
      cfg = merge_overrides(base, toml)
      cfg = merge_overrides(cfg,  env)
      cfg = merge_overrides(cfg,  capi)
   end function resolve_diagnostics_config

   !> Built-in default config for a standalone CLI run.
   pure function default_cli_config() result(cfg)
      type(diagnostics_config_t) :: cfg
      cfg%level      = LOGLEVEL_INFO
      cfg%to_stdout  = .true.
      cfg%to_stderr  = .true.
      cfg%timestamps = .false.
      cfg%log_file   = 'swap_swap.log'
   end function default_cli_config

   !> Built-in default config when embedded (BMI/C-API/XMI): never write the
   !! host's stdout; no log file unless the caller asks (env/C-API).
   pure function default_embedded_config() result(cfg)
      type(diagnostics_config_t) :: cfg
      cfg%level      = LOGLEVEL_INFO
      cfg%to_stdout  = .false.
      cfg%to_stderr  = .true.
      cfg%timestamps = .false.
      ! log_file deliberately left unallocated
   end function default_embedded_config

   !> Read logging overrides from the environment:
   !!   SWAP_LOG_LEVEL   = DEBUG|INFO|WARN|ERROR|NONE
   !!   SWAP_LOG_FILE    = <path>
   !!   SWAP_LOG_STDOUT  = 0|1 (or true/false)
   !!   SWAP_LOG_STDERR  = 0|1 (or true/false)
   !! Impure (reads the environment).
   subroutine read_env_overrides(ov)
      type(diag_overrides_t), intent(out) :: ov
      character(len=256) :: buf
      integer :: ln, st
      call get_environment_variable('SWAP_LOG_LEVEL', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_level = .true.
         ov%level     = log_level_from_name(buf(:ln))
      end if
      call get_environment_variable('SWAP_LOG_FILE', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_file = .true.
         ov%log_file = buf(:ln)
      end if
      call get_environment_variable('SWAP_LOG_STDOUT', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_stdout = .true.
         ov%to_stdout  = (buf(1:1) == '1' .or. buf(1:1) == 't' .or. buf(1:1) == 'T')
      end if
      call get_environment_variable('SWAP_LOG_STDERR', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_stderr = .true.
         ov%to_stderr  = (buf(1:1) == '1' .or. buf(1:1) == 't' .or. buf(1:1) == 'T')
      end if
   end subroutine read_env_overrides

   !> ASCII upper-case helper (pure, no locale).
   pure function upcase(s) result(u)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: u
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) then
            u(i:i) = achar(c - 32)
         else
            u(i:i) = s(i:i)
         end if
      end do
   end function upcase

end module diagnostics_mod
