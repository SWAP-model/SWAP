!> @file config_source.f90
!! A companion-content provider that abstracts *where* SWAP's companion
!! files (drainage `.dra.toml`, crop `.crp.toml`, CSV data tables) come
!! from. Standalone runs resolve names against a base directory on disk;
!! Python-driven runs resolve names against in-memory blobs registered by
!! the caller. Both yield identical text, so everything downstream of the
!! config loader is oblivious to the source — preserving byte-identical
!! behavior between the disk and in-memory paths.
!!
!! Threaded as an argument through the config-load chain (no module global).
module config_source_mod
   implicit none
   private

   public :: config_source_t
   public :: config_source_disk
   public :: config_source_memory

   type :: named_blob_t
      character(len=:), allocatable :: name
      character(len=:), allocatable :: content
   end type named_blob_t

   type :: config_source_t
      logical                       :: in_memory = .false.
      character(len=:), allocatable :: base_dir
      type(named_blob_t), allocatable :: blobs(:)
   contains
      procedure :: add_blob => config_source_add_blob
      procedure :: get_text => config_source_get_text
   end type config_source_t

contains

   !> Construct a disk-backed source rooted at `base_dir`.
   function config_source_disk(base_dir) result(src)
      character(len=*), intent(in) :: base_dir
      type(config_source_t)        :: src
      src%in_memory = .false.
      src%base_dir  = base_dir
      allocate (src%blobs(0))
   end function config_source_disk

   !> Construct an in-memory source; register companions with add_blob.
   function config_source_memory() result(src)
      type(config_source_t) :: src
      src%in_memory = .true.
      allocate (src%blobs(0))
   end function config_source_memory

   subroutine config_source_add_blob(self, name, content)
      class(config_source_t), intent(inout) :: self
      character(len=*),       intent(in)    :: name
      character(len=*),       intent(in)    :: content
      type(named_blob_t), allocatable :: grown(:)
      integer :: n

      if (.not. allocated(self%blobs)) allocate (self%blobs(0))
      n = size(self%blobs)
      allocate (grown(n + 1))
      grown(1:n) = self%blobs
      grown(n + 1)%name    = trim(name)
      grown(n + 1)%content = content
      call move_alloc(grown, self%blobs)
   end subroutine config_source_add_blob

   subroutine config_source_get_text(self, name, text, found)
      class(config_source_t),        intent(in)  :: self
      character(len=*),              intent(in)  :: name
      character(len=:), allocatable, intent(out) :: text
      logical,                       intent(out) :: found
      integer :: i

      found = .false.
      text  = ''

      if (self%in_memory) then
         if (allocated(self%blobs)) then
            do i = 1, size(self%blobs)
               if (self%blobs(i)%name == trim(name)) then
                  text  = self%blobs(i)%content
                  found = .true.
                  return
               end if
            end do
         end if
      else
         if (allocated(self%base_dir)) then
            call read_file_to_string(self%base_dir // trim(name), text, found)
         else
            call read_file_to_string(trim(name), text, found)
         end if
      end if
   end subroutine config_source_get_text

   !> Read an entire file into an allocatable string. found=.false. if the
   !! file cannot be opened. Stream access preserves the file's bytes
   !! (incl. newlines) so a TOML parser sees exactly the on-disk content.
   subroutine read_file_to_string(path, text, found)
      character(len=*),              intent(in)  :: path
      character(len=:), allocatable, intent(out) :: text
      logical,                       intent(out) :: found
      integer :: u, ios, fsize

      found = .false.
      open (newunit=u, file=path, access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
      if (ios /= 0) then
         text = ''
         return
      end if
      inquire (unit=u, size=fsize)
      if (fsize > 0) then
         allocate (character(len=fsize) :: text)
         read (u, iostat=ios) text
         if (ios /= 0) then
            close (u)
            text = ''
            return
         end if
      else
         text = ''
      end if
      close (u)
      found = .true.
   end subroutine read_file_to_string

end module config_source_mod
