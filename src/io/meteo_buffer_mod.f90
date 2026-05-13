!> @file meteo_buffer_mod.f90
!! SS-BMI2: external meteo buffer storage. When mode is set to
!! METEO_MODE_EXTERNAL_BUFFER (via swap_attach_meteo_buffer in the
!! cffi surface — Task 17), readmeteo skips its CSV-reading path
!! and reads from the externally-supplied buffer instead.
!!
!! Canonical column order:
!!   1: date_offset (days since simulation start, integer)
!!   2: rain        (mm/d)
!!   3: tmin        (deg C)
!!   4: tmax        (deg C)
!!   5: et_ref      (mm/d, -99.9 = compute internally)
!!   6: radiation   (kJ/m2/d)
!!   7: vapor       (kPa)
!!   8: wind        (m/s)
module meteo_buffer_mod
   use iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
   implicit none
   private

   public :: METEO_MODE_PATH, METEO_MODE_EXTERNAL_BUFFER
   public :: get_meteo_mode, set_meteo_mode_path
   public :: attach_external_meteo_buffer
   public :: get_external_meteo_value
   public :: get_external_meteo_n_days, get_external_meteo_n_cols

   integer, parameter :: METEO_MODE_PATH             = 0
   integer, parameter :: METEO_MODE_EXTERNAL_BUFFER  = 1

   integer,                save :: current_mode = METEO_MODE_PATH
   real(c_double), pointer      :: ext_buffer(:,:) => null()
   integer,                save :: ext_n_days = 0
   integer,                save :: ext_n_cols = 0

contains

   function get_meteo_mode() result(mode)
      integer :: mode
      mode = current_mode
   end function get_meteo_mode

   subroutine set_meteo_mode_path()
      current_mode = METEO_MODE_PATH
      nullify(ext_buffer)
      ext_n_days = 0
      ext_n_cols = 0
   end subroutine set_meteo_mode_path

   subroutine attach_external_meteo_buffer(ptr, n_days, n_cols)
      type(c_ptr),    value, intent(in) :: ptr
      integer(c_int), value, intent(in) :: n_days, n_cols
      call c_f_pointer(ptr, ext_buffer, [n_cols, n_days])
      ext_n_days = n_days
      ext_n_cols = n_cols
      current_mode = METEO_MODE_EXTERNAL_BUFFER
   end subroutine attach_external_meteo_buffer

   !> Read a single value from the attached buffer.
   !! @param day  1-based day index within the buffer
   !! @param col  1-based column index (canonical order — see module doc)
   function get_external_meteo_value(day, col) result(val)
      integer, intent(in) :: day, col
      real(c_double)      :: val
      val = ext_buffer(col, day)
   end function get_external_meteo_value

   function get_external_meteo_n_days() result(n)
      integer :: n
      n = ext_n_days
   end function get_external_meteo_n_days

   function get_external_meteo_n_cols() result(n)
      integer :: n
      n = ext_n_cols
   end function get_external_meteo_n_cols

end module meteo_buffer_mod
