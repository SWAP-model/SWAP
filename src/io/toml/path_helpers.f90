!> Path-string helpers for cross-file TOML references.
!!
!! Pure string manipulation — no filesystem calls. The loader uses
!! `directory_of` to extract a base directory from a TOML file path,
!! then `resolve_relative_path` to resolve cross-file references against
!! that base.
module path_helpers_mod
   implicit none
   private

   public :: directory_of
   public :: resolve_relative_path

contains

   !> Return the directory portion of `path`, including the trailing slash.
   !! If `path` has no slash, returns "./".
   function directory_of(path) result(dir)
      character(len=*),              intent(in)  :: path
      character(len=:), allocatable              :: dir
      integer :: i
      i = index(path, '/', back=.true.)
      if (i == 0) then
         dir = './'
      else
         dir = path(1:i)
      end if
   end function directory_of

   !> Resolve `rel` against `base`. If `rel` starts with `/` it is treated
   !! as absolute and returned verbatim. Otherwise returns `base // rel`.
   !! `base` is expected to end with a slash (use `directory_of`).
   function resolve_relative_path(base, rel) result(abs)
      character(len=*),              intent(in)  :: base, rel
      character(len=:), allocatable              :: abs
      if (len_trim(rel) > 0 .and. rel(1:1) == '/') then
         abs = trim(rel)
      else
         abs = trim(base) // trim(rel)
      end if
   end function resolve_relative_path

end module path_helpers_mod
