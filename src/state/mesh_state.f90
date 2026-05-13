!> @file mesh_state.f90
!! SS-GR-BH: typed mesh / vertical-discretization state record.
!! Replaces the bare globals `numnod`, `layer(:)`, `dz(:)`, `z(:)`,
!! `disnod(:)`, `ztopcp(:)`, `zbotcp(:)` from `variables.f90`.
!! Populated once by `mesh_init` from `config_to_variables.f90`;
!! treated as read-only thereafter.
module mesh_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: mesh_state_t

   type :: mesh_state_t
      integer :: numnod = 0
      integer,      allocatable :: layer(:)
      real(real64), allocatable :: dz(:)
      real(real64), allocatable :: z(:)
      real(real64), allocatable :: disnod(:)
      real(real64), allocatable :: ztopcp(:)
      real(real64), allocatable :: zbotcp(:)
   contains
      procedure :: init => mesh_init
   end type mesh_state_t

contains

   subroutine mesh_init(self, numnod_in, dz_in, z_in, disnod_in, &
                        ztopcp_in, zbotcp_in, layer_in)
      class(mesh_state_t), intent(inout) :: self
      integer,             intent(in)    :: numnod_in
      real(real64),        intent(in)    :: dz_in(:), z_in(:), disnod_in(:), &
                                            ztopcp_in(:), zbotcp_in(:)
      integer,             intent(in)    :: layer_in(:)

      self%numnod = numnod_in

      if (allocated(self%dz))     deallocate(self%dz)
      if (allocated(self%z))      deallocate(self%z)
      if (allocated(self%disnod)) deallocate(self%disnod)
      if (allocated(self%ztopcp)) deallocate(self%ztopcp)
      if (allocated(self%zbotcp)) deallocate(self%zbotcp)
      if (allocated(self%layer))  deallocate(self%layer)

      allocate(self%dz(numnod_in))
      allocate(self%z(numnod_in))
      allocate(self%disnod(numnod_in+1))
      allocate(self%ztopcp(numnod_in))
      allocate(self%zbotcp(numnod_in))
      allocate(self%layer(numnod_in))

      self%dz(:)     = dz_in(1:numnod_in)
      self%z(:)      = z_in(1:numnod_in)
      self%disnod(:) = disnod_in(1:numnod_in+1)
      self%ztopcp(:) = ztopcp_in(1:numnod_in)
      self%zbotcp(:) = zbotcp_in(1:numnod_in)
      self%layer(:)  = layer_in(1:numnod_in)
   end subroutine mesh_init

end module mesh_state_mod
