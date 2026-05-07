!> Reusable TOML field-reading primitives with error accumulation.
!!
!! Every helper wraps a tomlf get_value call, handles null-pointer
!! cases, and appends a typed error (ERR_PARSE_*) to the collection
!! on failure. Callers never check stat/associated themselves.
module toml_field_helpers_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value
   use error_mod, only: error_collection_t,            &
                        ERR_PARSE_TYPE_MISMATCH,       &
                        ERR_PARSE_MISSING_REQUIRED
   implicit none
   private

   public :: get_required_int
   public :: get_required_real
   public :: get_required_string
   public :: get_optional_int_with_default
   public :: get_optional_real_with_default
   public :: get_optional_logical_with_default
   public :: get_optional_string_with_default
   public :: get_table
   public :: get_array_of_tables
   public :: parse_date_to_days1900

contains

   subroutine get_required_int(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      integer,                   intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      if (.not. associated(tab)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "parent table missing for " // key, context)
         out = 0
         return
      end if
      call get_value(tab, key, out, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected integer at " // key, context)
         out = 0
      end if
   end subroutine get_required_int

   subroutine get_required_real(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      real(real64),              intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      if (.not. associated(tab)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "parent table missing for " // key, context)
         out = 0.0_real64
         return
      end if
      call get_value(tab, key, out, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected real at " // key, context)
         out = 0.0_real64
      end if
   end subroutine get_required_real

   subroutine get_required_string(tab, key, out, context, errors)
      type(toml_table), pointer,     intent(in)    :: tab
      character(len=*),              intent(in)    :: key
      character(len=:), allocatable, intent(out)   :: out
      character(len=*),              intent(in)    :: context
      type(error_collection_t),      intent(inout) :: errors
      integer :: stat
      character(len=:), allocatable :: tmp
      if (.not. associated(tab)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "parent table missing for " // key, context)
         out = ""
         return
      end if
      call get_value(tab, key, tmp, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected string at " // key, context)
         out = ""
         return
      end if
      out = tmp
   end subroutine get_required_string

   subroutine get_optional_int_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      integer,                   intent(out)   :: out
      integer,                   intent(in)    :: default
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, out, default=default, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected integer at " // key, context)
         out = default
      end if
   end subroutine get_optional_int_with_default

   subroutine get_optional_real_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      real(real64),              intent(out)   :: out
      real(real64),              intent(in)    :: default
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, out, default=default, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected real at " // key, context)
         out = default
      end if
   end subroutine get_optional_real_with_default

   subroutine get_optional_logical_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      logical,                   intent(out)   :: out
      logical,                   intent(in)    :: default
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, out, default=default, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected logical at " // key, context)
         out = default
      end if
   end subroutine get_optional_logical_with_default

   subroutine get_optional_string_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer,     intent(in)    :: tab
      character(len=*),              intent(in)    :: key
      character(len=:), allocatable, intent(out)   :: out
      character(len=*),              intent(in)    :: default
      character(len=*),              intent(in)    :: context
      type(error_collection_t),      intent(inout) :: errors
      integer :: stat
      character(len=:), allocatable :: tmp
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, tmp, stat=stat)
      if (stat /= 0 .or. .not. allocated(tmp)) then
         out = default
         return
      end if
      out = tmp
   end subroutine get_optional_string_with_default

   !> Look up a sub-table. Returns null on absence (caller decides whether
   !! absence is an error).
   subroutine get_table(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      type(toml_table), pointer, intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out => null()
      if (.not. associated(tab)) return
      call get_value(tab, key, out, requested=.false., stat=stat)
   end subroutine get_table

   !> Fetch an array-of-tables under key. Null on absence.
   subroutine get_array_of_tables(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      type(toml_array), pointer, intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out => null()
      if (.not. associated(tab)) return
      call get_value(tab, key, out, requested=.false., stat=stat)
   end subroutine get_array_of_tables

   !> Convert a TOML datetime to days-since-1900 (the legacy time axis).
   !! 1900-01-01 is JD 2415021 (but we use 2415020 for 1-based indexing).
   function parse_date_to_days1900(dtv) result(days)
      type(toml_datetime), intent(in) :: dtv
      real(real64) :: days
      integer :: y, m, d, jd, jd1900
      y = dtv%date%year
      m = dtv%date%month
      d = dtv%date%day
      jd     = julian_day(y, m, d)
      jd1900 = 2415020
      days   = real(jd - jd1900, kind=real64)
   end function parse_date_to_days1900

   pure function julian_day(y, m, d) result(jd)
      integer, intent(in) :: y, m, d
      integer :: jd, a, yy, mm
      a  = (14 - m) / 12
      yy = y + 4800 - a
      mm = m + 12 * a - 3
      jd = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045
   end function julian_day

end module toml_field_helpers_mod
