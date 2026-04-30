!> Error handling for SWAP.
!!
!! Defines a typed error payload (`error_t`) and an accumulating
!! collection (`error_collection_t`). Every fallible procedure in the
!! new infrastructure takes `errors` as `intent(inout)` and appends
!! failures without aborting. A single abort point (`abort_if_fatal`)
!! turns accumulated fatal errors into `error stop`.
module error_mod
   use iso_fortran_env, only: real64, error_unit
   use swap_log, only: log_error
   implicit none
   private

   ! Error code parameters. Stable values; tests assert on these.
   integer, parameter, public :: ERR_NONE                     = 0
   integer, parameter, public :: ERR_IO_READ_FAILED           = 100
   integer, parameter, public :: ERR_IO_WRITE_FAILED          = 101
   integer, parameter, public :: ERR_IO_OPEN_FAILED           = 102
   integer, parameter, public :: ERR_PARSE_MALFORMED_TOML     = 200
   integer, parameter, public :: ERR_PARSE_TYPE_MISMATCH      = 201
   integer, parameter, public :: ERR_PARSE_MISSING_REQUIRED   = 202
   integer, parameter, public :: ERR_PARSE_MISSING_HEADER     = 203
   integer, parameter, public :: ERR_PARSE_HEADER_MISMATCH    = 204
   integer, parameter, public :: ERR_PARSE_ROW_SHAPE          = 205
   integer, parameter, public :: ERR_VALIDATION_OUT_OF_RANGE  = 300
   integer, parameter, public :: ERR_VALIDATION_ENUM          = 301
   integer, parameter, public :: ERR_VALIDATION_CROSS_FIELD   = 302
   integer, parameter, public :: ERR_VALIDATION_CROSS_SECTION = 303
   integer, parameter, public :: ERR_VALIDATION_REQUIRED      = 304
   integer, parameter, public :: ERR_FINALIZE_DERIVATION      = 400
   integer, parameter, public :: ERR_ADAPTER_UNSUPPORTED      = 500
   integer, parameter, public :: ERR_LEGACY_FATAL             = 999
   integer, parameter, public :: ERR_DEPRECATED_KEY            = 600

   type, public :: error_t
      integer                       :: code    = ERR_NONE
      character(len=:), allocatable :: message
      character(len=:), allocatable :: context
      logical                       :: is_fatal = .false.
   end type error_t

   type, public :: error_collection_t
      type(error_t), allocatable :: items(:)
   contains
      procedure :: append         => error_collection_append
      procedure :: has_errors     => error_collection_has_errors
      procedure :: has_fatals     => error_collection_has_fatals
      procedure :: count          => error_collection_count
      procedure :: summary        => error_collection_summary
      procedure :: abort_if_fatal => error_collection_abort_if_fatal
      procedure :: clear          => error_collection_clear
   end type error_collection_t

   !> Module-level singleton — drop-in for legacy `fatalerr` calls deep
   !! in physics paths where threading a per-call `errors` argument is
   !! impractical. Top-level test setup may call `clear()` between runs.
   type(error_collection_t), public, save :: global_errors

   public :: fatalerr_collected
   public :: warn_deprecated_key

contains

   !> Append a non-fatal deprecation warning for a TOML key that was
   !! retired per an ADR (e.g. ADR 0009 retired the legacy non-CSV
   !! output switches). Used by new TOML readers when they encounter
   !! a key that has no schema slot but appears in user input.
   !! The warning is NOT fatal — the reader continues and ignores
   !! the key. Use abort_if_fatal at the pipeline edge if a downstream
   !! check decides the deprecation is now fatal.
   subroutine warn_deprecated_key(routine, key, errors)
      character(len=*),         intent(in)    :: routine
      character(len=*),         intent(in)    :: key
      type(error_collection_t), intent(inout) :: errors
      character(len=256) :: msg
      write(msg, '("deprecated key ''", A, "'' is ignored; see ADR 0009")') trim(key)
      call errors%append(ERR_DEPRECATED_KEY, trim(msg), routine, is_fatal=.false.)
   end subroutine warn_deprecated_key

   !> Drop-in replacement for legacy `call fatalerr(routine, msg)`.
   !! Appends a fatal entry to `global_errors` then aborts. Used by
   !! physics-path code where threading an explicit `errors` argument
   !! through every caller would be invasive. I/O and configuration
   !! code should still use the threaded-`errors` pattern instead.
   subroutine fatalerr_collected(routine, message)
      character(len=*), intent(in) :: routine
      character(len=*), intent(in) :: message
      call global_errors%append(ERR_LEGACY_FATAL, message, routine)
      call global_errors%abort_if_fatal()
   end subroutine fatalerr_collected

   !> Append an error to the collection and auto-log through swap_log.
   subroutine error_collection_append(self, code, message, context, is_fatal)
      class(error_collection_t), intent(inout) :: self
      integer,                   intent(in)    :: code
      character(len=*),          intent(in)    :: message
      character(len=*),          intent(in)    :: context
      logical, optional,         intent(in)    :: is_fatal

      type(error_t), allocatable :: tmp(:)
      type(error_t)              :: item
      integer                    :: n

      item%code     = code
      item%message  = message
      item%context  = context
      if (present(is_fatal)) then
         item%is_fatal = is_fatal
      else
         item%is_fatal = .true.
      end if

      if (.not. allocated(self%items)) then
         allocate(self%items(1))
         self%items(1) = item
      else
         n = size(self%items)
         allocate(tmp(n + 1))
         tmp(1:n)   = self%items
         tmp(n + 1) = item
         call move_alloc(from=tmp, to=self%items)
      end if

      call log_error(context, message)
   end subroutine error_collection_append

   !> .true. if any items have been appended.
   pure function error_collection_has_errors(self) result(yes)
      class(error_collection_t), intent(in) :: self
      logical :: yes
      yes = allocated(self%items) .and. size_safe(self) > 0
   end function error_collection_has_errors

   !> .true. if any appended item has is_fatal=.true.
   pure function error_collection_has_fatals(self) result(yes)
      class(error_collection_t), intent(in) :: self
      logical :: yes
      integer :: i
      yes = .false.
      if (.not. allocated(self%items)) return
      do i = 1, size(self%items)
         if (self%items(i)%is_fatal) then
            yes = .true.
            return
         end if
      end do
   end function error_collection_has_fatals

   !> Number of errors appended so far.
   pure function error_collection_count(self) result(n)
      class(error_collection_t), intent(in) :: self
      integer :: n
      if (allocated(self%items)) then
         n = size(self%items)
      else
         n = 0
      end if
   end function error_collection_count

   !> Multi-line human-readable report of all collected errors.
   function error_collection_summary(self) result(text)
      class(error_collection_t), intent(in) :: self
      character(len=:), allocatable :: text
      character(len=32)             :: code_str
      integer                       :: i

      if (.not. allocated(self%items) .or. size(self%items) == 0) then
         text = "No errors."
         return
      end if

      text = ""
      do i = 1, size(self%items)
         write(code_str, '(I0)') self%items(i)%code
         text = text // "[" // trim(adjustl(code_str)) // "]"
         if (self%items(i)%is_fatal) then
            text = text // " FATAL "
         else
            text = text // " WARN  "
         end if
         text = text // self%items(i)%context // ": " // self%items(i)%message // new_line('a')
      end do
   end function error_collection_summary

   !> Write summary to stderr and `error stop` if any fatal errors exist.
   !! Returns normally if no fatals.
   subroutine error_collection_abort_if_fatal(self)
      class(error_collection_t), intent(in) :: self
      if (.not. self%has_fatals()) return
      write(error_unit, '(A)') self%summary()
      error stop "fatal error(s) in swap input pipeline"
   end subroutine error_collection_abort_if_fatal

   !> Reset the collection. Use between tests or logical phases.
   subroutine error_collection_clear(self)
      class(error_collection_t), intent(inout) :: self
      if (allocated(self%items)) deallocate(self%items)
   end subroutine error_collection_clear

   !> Helper used by has_errors (keeps the function pure).
   pure function size_safe(self) result(n)
      class(error_collection_t), intent(in) :: self
      integer :: n
      if (allocated(self%items)) then
         n = size(self%items)
      else
         n = 0
      end if
   end function size_safe

end module error_mod
