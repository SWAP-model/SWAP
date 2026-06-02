!> @file results_state.f90
!! Per-instance in-memory record of SWAP's scalar output. Holds only the
!! user's selected (inlist) columns. Rows are appended via add_row and
!! committed to a growing in-memory store by flush every `flush_every` rows;
!! the chunk/flush seam is the swap-in point for a future binary disk sink.
!! Scalar output only (depth-time profile output stays on its streaming path).
module results_state_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: results_state_t

   type :: results_state_t
      integer :: ncols       = 0
      integer :: nrows       = 0          !! committed rows
      integer :: flush_every = 512        !! rows per chunk before a flush
      character(len=24), allocatable :: col_names(:)
      character(len=12), allocatable :: col_units(:)
      real(real64),      allocatable :: times(:)        !! (nrows)
      real(real64),      allocatable :: values(:,:)     !! (nrows, ncols)
      ! chunk buffer (filled by add_row, drained by flush)
      integer :: chunk_n = 0
      real(real64), allocatable :: chunk_t(:)           !! (flush_every)
      real(real64), allocatable :: chunk_v(:,:)         !! (flush_every, ncols)
   contains
      procedure :: init     => results_state_init
      procedure :: add_row  => results_state_add_row
      procedure :: flush    => results_state_flush
      procedure :: finalize => results_state_finalize
   end type results_state_t

contains

   subroutine results_state_init(self, names, units, flush_every)
      class(results_state_t), intent(inout) :: self
      character(len=*),       intent(in)    :: names(:)
      character(len=*),       intent(in)    :: units(:)
      integer, optional,      intent(in)    :: flush_every

      self%ncols = size(names)
      self%nrows = 0
      self%chunk_n = 0
      if (present(flush_every)) self%flush_every = max(1, flush_every)
      self%col_names = names
      self%col_units = units
      allocate(self%times(0))
      allocate(self%values(0, self%ncols))
      allocate(self%chunk_t(self%flush_every))
      allocate(self%chunk_v(self%flush_every, self%ncols))
   end subroutine results_state_init

   subroutine results_state_add_row(self, time, vals)
      class(results_state_t), intent(inout) :: self
      real(real64),           intent(in)    :: time
      real(real64),           intent(in)    :: vals(:)
      ! Caller invariant: exactly one value per configured column. A mismatch
      ! means the column set and the emitted row drifted apart — a programming
      ! bug we surface loudly rather than silently truncate/zero-pad.
      if (size(vals) /= self%ncols) &
         error stop 'results_state_t%add_row: size(vals) /= ncols'
      self%chunk_n = self%chunk_n + 1
      self%chunk_t(self%chunk_n)               = time
      self%chunk_v(self%chunk_n, 1:self%ncols) = vals
      if (self%chunk_n == self%flush_every) call self%flush()
   end subroutine results_state_add_row

   !> Drain the chunk buffer into the committed store. For the in-memory sink
   !! this grows `times`/`values`; a future disk sink would write the chunk to
   !! disk here and leave the in-memory store empty.
   subroutine results_state_flush(self)
      class(results_state_t), intent(inout) :: self
      real(real64), allocatable :: t2(:), v2(:,:)
      integer :: old, add
      if (self%chunk_n == 0) return
      old = self%nrows
      add = self%chunk_n
      allocate(t2(old + add))
      allocate(v2(old + add, self%ncols))
      if (old > 0) then
         t2(1:old)      = self%times
         v2(1:old, :)   = self%values
      end if
      t2(old+1:old+add)    = self%chunk_t(1:add)
      v2(old+1:old+add, :) = self%chunk_v(1:add, :)
      call move_alloc(t2, self%times)
      call move_alloc(v2, self%values)
      self%nrows   = old + add
      self%chunk_n = 0
   end subroutine results_state_flush

   subroutine results_state_finalize(self)
      class(results_state_t), intent(inout) :: self
      call self%flush()
   end subroutine results_state_finalize

end module results_state_mod
