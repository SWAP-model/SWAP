!> Generic CSV writer for SWAP result files — symmetric counterpart of
!! csv_reader_mod. A csv_writer_t owns one open unit. The number formatter
!! (fmt_real) preserves the legacy E/F selection from swap_csv_output.f90's
!! what_form so byte output is stable.
!!
!! Formatter spec (lifted verbatim from what_form, DO NOT re-derive):
!!   num_d = 5  (decimals)
!!   num_w = 12 (= num_d + 7, total E-field width)
!!   expo  = 4  => t1 = 1.0d-4, t2 = 1.0d4
!!   |x| < 1.0D-10          => form_rea_F0 = '(F0.0)'   (zero decimals)
!!   |x| < t1 or |x| > t2  => form_rea_E1 = '(1P,E12.5)'
!!   otherwise              => form_rea_F1 = '(F0.5)'
module csv_writer_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t, ERR_IO_OPEN_FAILED
   use file_io_mod,     only: file_open
   implicit none
   private

   public :: csv_writer_t
   public :: fmt_real_for_test   ! thin alias exported only for unit tests

   type :: csv_writer_t
      integer :: unit  = -1
      integer :: ncols = 0
   contains
      procedure :: open   => csv_writer_open
      procedure :: meta   => csv_writer_meta
      procedure :: header => csv_writer_header
      procedure :: row    => csv_writer_row
      procedure :: flush  => csv_writer_flush
      procedure :: close  => csv_writer_close
   end type csv_writer_t

contains

   subroutine csv_writer_open(self, path, errors)
      class(csv_writer_t),      intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors
      integer :: ios
      call file_open(self%unit, path, 'replace', 'write', iostat=ios)
      if (ios /= 0) then
         call errors%append(ERR_IO_OPEN_FAILED, &
            "cannot open CSV for write: " // trim(path), 'csv_writer')
         self%unit = -1
      end if
   end subroutine csv_writer_open

   !> Write '*'-prefixed metadata comment lines (Project, File content, ...).
   subroutine csv_writer_meta(self, lines)
      class(csv_writer_t), intent(inout) :: self
      character(len=*),    intent(in)    :: lines(:)
      integer :: i
      if (self%unit == -1) return
      do i = 1, size(lines)
         write(self%unit, '(A)') '* ' // trim(lines(i))
      end do
   end subroutine csv_writer_meta

   !> Write the column-name row and (optionally) a unit row. Sets ncols.
   subroutine csv_writer_header(self, names, units)
      class(csv_writer_t), intent(inout) :: self
      character(len=*),    intent(in)    :: names(:)
      character(len=*),    intent(in), optional :: units(:)
      if (self%unit == -1) return
      self%ncols = size(names)
      write(self%unit, '(A)') join(names)
      if (present(units)) write(self%unit, '(A)') join(units)
   end subroutine csv_writer_header

   !> Write one data row: optional leading string (datetime), then values.
   subroutine csv_writer_row(self, values, leading)
      class(csv_writer_t), intent(inout) :: self
      real(real64),        intent(in)    :: values(:)
      character(len=*),    intent(in), optional :: leading
      character(len=:), allocatable :: line
      integer :: j
      if (self%unit == -1) return
      line = ''
      if (present(leading)) line = trim(leading) // ','
      do j = 1, size(values)
         line = line // trim(adjustl(fmt_real(values(j))))
         if (j < size(values)) line = line // ','
      end do
      write(self%unit, '(A)') line
   end subroutine csv_writer_row

   !> Drain the runtime I/O buffer to disk without closing. No-op if not open.
   subroutine csv_writer_flush(self)
      class(csv_writer_t), intent(inout) :: self
      if (self%unit /= -1) flush(self%unit)
   end subroutine csv_writer_flush

   subroutine csv_writer_close(self)
      class(csv_writer_t), intent(inout) :: self
      if (self%unit /= -1) close(self%unit)
      self%unit = -1
   end subroutine csv_writer_close

   ! ---- private helpers --------------------------------------------------

   pure function join(fields) result(s)
      character(len=*), intent(in) :: fields(:)
      character(len=:), allocatable :: s
      integer :: i
      s = ''
      do i = 1, size(fields)
         s = s // trim(adjustl(fields(i)))
         if (i < size(fields)) s = s // ','
      end do
   end function join

   !> Legacy E/F selector lifted verbatim from swap_csv_output.f90's what_form.
   !! Thresholds and format strings must NOT be re-derived — the regression
   !! tolerance (1e-2) depends on byte-stable formatting.
   !!
   !! what_form parameters (source of truth):
   !!   num_d = 5, num_w = 12, expo = 4
   !!   t1 = 1.0d0/(10.0d0**4) = 1.0d-4
   !!   t2 = 10.0d0**4          = 1.0d4
   !!   form_rea_F0 = '(F0.0,",")'   when |x| < 1.0D-10
   !!   form_rea_E1 = '(1P,E12.5,",")'  when |x|<t1 or |x|>t2
   !!   form_rea_F1 = '(F0.5,",")'   otherwise
   !! (comma suffix stripped here; caller joins with ',')
   function fmt_real(x) result(buf)
      real(real64), intent(in) :: x
      character(len=30) :: buf
      real(real64), parameter :: t1 = 1.0d0 / (10.0d0**4)
      real(real64), parameter :: t2 = 10.0d0**4
      if (abs(x) < 1.0d-10) then
         write(buf, '(F0.0)') x
      else if (abs(x) < t1 .or. abs(x) > t2) then
         write(buf, '(1P,E12.5)') x
      else
         write(buf, '(F0.5)') x
      end if
   end function fmt_real

   function fmt_real_for_test(x) result(buf)
      real(real64), intent(in) :: x
      character(len=30) :: buf
      buf = fmt_real(x)
   end function fmt_real_for_test

end module csv_writer_mod
