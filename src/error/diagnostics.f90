!> Diagnostics configuration layer (Phase A).
!!
!! Pure resolution of logging settings from layered sources with the
!! precedence  C-API > env > TOML > built-in default. This
!! phase owns only the *config* (level/file/routing) + a precedence merge;
!! the per-instance state-borne diagnostics record is a later phase.
module diagnostics_mod
   use swap_log, only: LOGLEVEL_DEBUG, LOGLEVEL_INFO, LOGLEVEL_WARN, &
                       LOGLEVEL_ERROR, LOGLEVEL_NONE, log_init, &
                       log_debug, log_info, log_warn, log_error
   use error_mod, only: error_collection_t, ERR_LEGACY_FATAL
   implicit none
   private

   public :: diagnostics_config_t, diag_overrides_t
   public :: diagnostics_t
   public :: log_level_from_name
   public :: merge_overrides, resolve_diagnostics_config
   public :: read_env_overrides, default_cli_config, default_embedded_config
   public :: init_logging_from_env, init_logging
   public :: read_logging_overrides_from_file, read_logging_overrides_from_text

   !> Per-instance diagnostics carried on swap_state_t (state%diag). Owns the
   !! instance's error accumulation + fatal flag + sim-time context; leveled
   !! emits route through the process-global swap_log sink, stamped with
   !! instance-id and sim-date so interleaved multi-instance logs are
   !! attributable.
   type, public :: diagnostics_t
      integer                  :: instance_id  = 0
      character(len=11)        :: sim_date     = ''
      integer                  :: daynr        = 0
      integer                  :: daycum       = 0
      type(error_collection_t) :: errors
      logical                  :: fatal_raised = .false.
   contains
      procedure :: debug       => diag_debug
      procedure :: info        => diag_info
      procedure :: warn        => diag_warn
      procedure :: error       => diag_error
      procedure :: fatal       => diag_fatal
      procedure :: aborted     => diag_aborted
      procedure :: set_simtime => diag_set_simtime
      procedure, private :: stamp => diag_stamp
   end type diagnostics_t

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

   !> Resolve a base config against the environment (env > base) and
   !! initialise the logger. Entry points call this with their default
   !! (CLI vs embedded). Higher-precedence C-API/TOML overrides are applied
   !! afterwards via their own setters.
   subroutine init_logging_from_env(base)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t)     :: none_ov
      call init_logging(base, none_ov)
   end subroutine init_logging_from_env

   !> Resolve base < toml < env and initialise the logger (one resolve, early).
   !! Entry points read their [logging] overrides and pass them here; the env
   !! still wins over TOML, and any higher-precedence C-API override is applied
   !! by its own setter afterwards.
   subroutine init_logging(base, toml_ov)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t),     intent(in) :: toml_ov
      type(diagnostics_config_t) :: dcfg
      type(diag_overrides_t)     :: none_ov, env_ov
      call read_env_overrides(env_ov)
      dcfg = resolve_diagnostics_config(base, toml_ov, env_ov, none_ov)
      call apply_logging_config(dcfg)
   end subroutine init_logging

   !> Hand a fully-resolved diagnostics config to the logger backend.
   subroutine apply_logging_config(dcfg)
      type(diagnostics_config_t), intent(in) :: dcfg
      if (allocated(dcfg%log_file)) then
         call log_init(log_level=dcfg%level, log_file=dcfg%log_file, &
                       to_stdout=dcfg%to_stdout, to_stderr=dcfg%to_stderr, &
                       timestamps=dcfg%timestamps)
      else
         call log_init(log_level=dcfg%level, to_stdout=dcfg%to_stdout, &
                       to_stderr=dcfg%to_stderr, timestamps=dcfg%timestamps)
      end if
   end subroutine apply_logging_config

   !> Extract [logging] overrides from a parsed TOML document.
   !! Recognised keys: level (string), file (string), to_stdout, to_stderr,
   !! timestamps (booleans). Absent [logging] table => empty overrides.
   subroutine extract_logging_overrides(doc, ov)
      use tomlf, only: toml_table, get_value
      type(toml_table), pointer, intent(in)  :: doc
      type(diag_overrides_t),    intent(out) :: ov
      type(toml_table), pointer     :: sec
      character(len=:), allocatable :: sval
      logical :: bval
      integer :: stat
      call get_value(doc, 'logging', sec, requested=.false.)
      if (.not. associated(sec)) return
      call get_value(sec, 'level', sval, stat=stat)
      if (stat == 0 .and. allocated(sval)) then
         if (len_trim(sval) > 0) then
            ov%has_level = .true.
            ov%level     = log_level_from_name(sval)
         end if
      end if
      call get_value(sec, 'file', sval, stat=stat)
      if (stat == 0 .and. allocated(sval)) then
         if (len_trim(sval) > 0) then
            ov%has_file  = .true.
            ov%log_file  = sval
         end if
      end if
      call get_value(sec, 'to_stdout', bval, stat=stat)
      if (stat == 0) then
         ov%has_stdout = .true.
         ov%to_stdout  = bval
      end if
      call get_value(sec, 'to_stderr', bval, stat=stat)
      if (stat == 0) then
         ov%has_stderr = .true.
         ov%to_stderr  = bval
      end if
      call get_value(sec, 'timestamps', bval, stat=stat)
      if (stat == 0) then
         ov%has_timestamps = .true.
         ov%timestamps     = bval
      end if
   end subroutine extract_logging_overrides

   !> Parse [logging] overrides from a TOML file path (missing/malformed => empty).
   subroutine read_logging_overrides_from_file(path, ov)
      use tomlf, only: toml_table, toml_error, toml_load
      character(len=*),       intent(in)  :: path
      type(diag_overrides_t), intent(out) :: ov
      type(toml_table), allocatable, target :: doc
      type(toml_table), pointer             :: doc_ptr
      type(toml_error), allocatable         :: terr
      call toml_load(doc, path, error=terr)
      if (allocated(terr)) return
      doc_ptr => doc
      call extract_logging_overrides(doc_ptr, ov)
   end subroutine read_logging_overrides_from_file

   !> Parse [logging] overrides from a TOML text buffer (in-memory C-API path).
   subroutine read_logging_overrides_from_text(text, ov)
      use tomlf, only: toml_table, toml_error, toml_loads
      character(len=*),       intent(in)  :: text
      type(diag_overrides_t), intent(out) :: ov
      type(toml_table), allocatable, target :: doc
      type(toml_table), pointer             :: doc_ptr
      type(toml_error), allocatable         :: terr
      call toml_loads(doc, text, error=terr)
      if (allocated(terr)) return
      doc_ptr => doc
      call extract_logging_overrides(doc_ptr, ov)
   end subroutine read_logging_overrides_from_text

   !> Prepend "[#id] date " to a message (id omitted when 0, date when empty)
   !! so single-instance CLI logs stay clean and ensemble logs are attributable.
   function diag_stamp(self, message) result(s)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: message
      character(len=:), allocatable    :: s
      character(len=16) :: idbuf
      s = message
      if (len_trim(self%sim_date) > 0) s = trim(self%sim_date)//' '//s
      if (self%instance_id /= 0) then
         write(idbuf,'("[#",I0,"] ")') self%instance_id
         s = trim(idbuf)//s
      end if
   end function diag_stamp

   subroutine diag_debug(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_debug(context, self%stamp(message))
   end subroutine diag_debug

   subroutine diag_info(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_info(context, self%stamp(message))
   end subroutine diag_info

   subroutine diag_warn(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_warn(context, self%stamp(message))
   end subroutine diag_warn

   subroutine diag_error(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_error(context, self%stamp(message))
   end subroutine diag_error

   !> Record a fatal into the per-instance collection (auto-logs via log_error)
   !! and set the sticky fatal flag. Does NOT abort — the step driver checks
   !! aborted() at substep boundaries (Phase D).
   subroutine diag_fatal(self, context, message)
      class(diagnostics_t), intent(inout) :: self
      character(len=*),     intent(in)    :: context, message
      call self%errors%append(ERR_LEGACY_FATAL, self%stamp(message), context)
      self%fatal_raised = .true.
   end subroutine diag_fatal

   pure logical function diag_aborted(self)
      class(diagnostics_t), intent(in) :: self
      diag_aborted = self%fatal_raised
   end function diag_aborted

   subroutine diag_set_simtime(self, date, daynr, daycum)
      class(diagnostics_t), intent(inout) :: self
      character(len=*),     intent(in)    :: date
      integer,              intent(in)    :: daynr, daycum
      self%sim_date = date
      self%daynr    = daynr
      self%daycum   = daycum
   end subroutine diag_set_simtime

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
