!> @file file_io.f90
!! Thin wrapper around Fortran intrinsic file I/O (open/close/inquire)
!! with swap_log integration. Replaces TTutil's getun/fopens/delfil.
!!
!! Conventions:
!! - file_open allocates a fresh unit via newunit=.
!! - status accepts native Fortran 'old' / 'new' / 'replace' /
!!   'unknown' / 'scratch'.
!! - action accepts native Fortran 'read' / 'write' / 'readwrite'.
!! - When iostat is provided, callers handle errors. When absent,
!!   open failure routes through fatalerr_collected (exit 1).
!! - Failures are logged at WARN level with path + iostat.
!!
!! See ADR 0023 for the umbrella context (TTutil retirement).
module file_io_mod
   use error_mod,  only: fatalerr_collected
   use swap_log,   only: log_warn, to_str
   implicit none
   private

   public :: file_open
   public :: file_delete
   public :: file_exists

contains

   !> Open `path` with the given Fortran status/action, allocating
   !! a new unit number into `unit`.
   !!
   !! On failure: when `iostat` is present, sets it nonzero and
   !! logs a WARN entry. When `iostat` is absent, logs and aborts
   !! via fatalerr_collected.
   subroutine file_open(unit, path, status, action, iostat)
      integer,           intent(out) :: unit
      character(len=*),  intent(in)  :: path
      character(len=*),  intent(in)  :: status
      character(len=*),  intent(in)  :: action
      integer, optional, intent(out) :: iostat
      integer :: ios

      open(newunit=unit, file=trim(path), status=trim(status), &
           action=trim(action), iostat=ios)

      if (ios /= 0) then
         call log_warn('file_io', "open failed: '" // trim(path) // &
              "' status=" // trim(status) // " action=" // trim(action) // &
              " iostat=" // trim(to_str(ios)))
         if (.not. present(iostat)) then
            call fatalerr_collected('file_io', &
                 "cannot open '" // trim(path) // &
                 "' (status=" // trim(status) // &
                 ", action=" // trim(action) // ")")
         end if
      end if

      if (present(iostat)) iostat = ios
   end subroutine file_open


   !> Delete `path` if it exists. Idempotent — no-op (and no error)
   !! if the file is missing.
   subroutine file_delete(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      logical :: exists

      inquire(file=trim(path), exist=exists)
      if (.not. exists) return

      open(newunit=u, file=trim(path), status='old', iostat=ios)
      if (ios /= 0) then
         call log_warn('file_io', "delete: cannot open '" // trim(path) // &
              "' iostat=" // trim(to_str(ios)))
         return
      end if
      close(u, status='delete')
   end subroutine file_delete


   !> True iff `path` refers to an existing file. Side-effect-free.
   logical function file_exists(path) result(exists)
      character(len=*), intent(in) :: path
      inquire(file=trim(path), exist=exists)
   end function file_exists

end module file_io_mod
