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
   integer, parameter, public :: ERR_PARSE_MALFORMED_TOML     = 200
   integer, parameter, public :: ERR_PARSE_TYPE_MISMATCH      = 201
   integer, parameter, public :: ERR_PARSE_MISSING_REQUIRED   = 202
   integer, parameter, public :: ERR_VALIDATION_OUT_OF_RANGE  = 300
   integer, parameter, public :: ERR_VALIDATION_ENUM          = 301
   integer, parameter, public :: ERR_VALIDATION_CROSS_FIELD   = 302
   integer, parameter, public :: ERR_VALIDATION_CROSS_SECTION = 303
   integer, parameter, public :: ERR_FINALIZE_DERIVATION      = 400
   integer, parameter, public :: ERR_ADAPTER_UNSUPPORTED      = 500

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

contains

   subroutine error_collection_append(self, code, message, context, is_fatal)
      class(error_collection_t), intent(inout) :: self
      integer,                   intent(in)    :: code
      character(len=*),          intent(in)    :: message
      character(len=*),          intent(in)    :: context
      logical, optional,         intent(in)    :: is_fatal
      ! stub - implementation in Task 3
   end subroutine error_collection_append

   pure function error_collection_has_errors(self) result(yes)
      class(error_collection_t), intent(in) :: self
      logical :: yes
      yes = .false.
   end function error_collection_has_errors

   pure function error_collection_has_fatals(self) result(yes)
      class(error_collection_t), intent(in) :: self
      logical :: yes
      yes = .false.
   end function error_collection_has_fatals

   pure function error_collection_count(self) result(n)
      class(error_collection_t), intent(in) :: self
      integer :: n
      n = 0
   end function error_collection_count

   function error_collection_summary(self) result(text)
      class(error_collection_t), intent(in) :: self
      character(len=:), allocatable :: text
      text = ""
   end function error_collection_summary

   subroutine error_collection_abort_if_fatal(self)
      class(error_collection_t), intent(in) :: self
      ! stub - implementation in Task 5
   end subroutine error_collection_abort_if_fatal

   subroutine error_collection_clear(self)
      class(error_collection_t), intent(inout) :: self
      ! stub - implementation in Task 4
   end subroutine error_collection_clear

end module error_mod
