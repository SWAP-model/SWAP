!> @file mesh_state.f90
!! SS-GR-BH: typed mesh / vertical-discretization state record.
!! Replaces the bare globals `numnod`, `layer(:)`, `dz(:)`, `z(:)`,
!! `disnod(:)`, `ztopcp(:)`, `zbotcp(:)` from `variables.f90`.
!! Populated once by `state%mesh%init(config%soil)` (formerly the free
!! subroutine `calcgrid`); treated as read-only thereafter.
module mesh_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use soil_config_mod,       only: soil_config_t
   use swap_array_dimensions, only: macp, maho
   use error_mod,             only: fatalerr_collected
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
      ! Linear-interpolation weights between neighbouring nodes (derived
      ! from dz/disnod). Filled by init; consumed by solute.
      real(real64), allocatable :: inpola(:)
      real(real64), allocatable :: inpolb(:)
      ! Per-layer indices computed by init (mesh-derived runtime state).
      integer :: numlay = 0                   !! number of physical soil layers
      integer, allocatable :: botcom(:)       !! bottom compartment of each soil layer
      integer, allocatable :: nod1lay(:)      !! first node of each soil layer
   contains
      procedure :: init => mesh_init
   end type mesh_state_t

contains

   !> Build the vertical discretization (formerly the free subroutine
   !> calcgrid). Computes grid geometry + layer/node maps from the soil
   !> config and fills this mesh record.
   subroutine mesh_init(self, soil_config)
      class(mesh_state_t), intent(inout) :: self
      type(soil_config_t), intent(in)    :: soil_config

      call validate_compartment_heights(soil_config)
      call allocate_grid_arrays(self)
      call place_nodes(self, soil_config)
      call compute_compartment_depths(self)
      call map_layers_to_nodes(self)
      call compute_interpolation_weights(self)
   end subroutine mesh_init

   !> Check correct input of number and height of soil compartments.
   subroutine validate_compartment_heights(soil_config)
      type(soil_config_t), intent(in) :: soil_config

      integer :: i, nsublay
      real(real64) :: hcomp_i
      character(len=200) :: messag
      character(len=11)  :: tmp

      ! sublay-derived counts come from the config arrays.
      nsublay = 0
      if (allocated(soil_config%sublay)) nsublay = size(soil_config%sublay)

      ! Check correct input of number and height of soil compartments.
      do i = 1, nsublay
         hcomp_i = soil_config%hsublay(i) / dble(soil_config%ncomp(i))
         if (abs(soil_config%ncomp(i)*hcomp_i - soil_config%hsublay(i)) .gt. 1.d-5) then
            write(tmp, '(i11)') i
            tmp = adjustl(tmp)
            messag = 'At the soil water section, part 4, at layer '//   &
         &           trim(tmp)//' the height of this soil layer hsublay corresponds'// &
         &           ' not to the product of height and number of compartments'
            call fatalerr_collected('mesh_init', messag)
         end if
      end do
   end subroutine validate_compartment_heights

   !> Allocate mesh arrays sized for macp (upper bound; trimmed below).
   subroutine allocate_grid_arrays(self)
      type(mesh_state_t), intent(inout) :: self

      ! The dealloc guard makes init idempotent / safe for re-init.
      if (allocated(self%dz))     deallocate(self%dz)
      if (allocated(self%z))      deallocate(self%z)
      if (allocated(self%disnod)) deallocate(self%disnod)
      if (allocated(self%ztopcp)) deallocate(self%ztopcp)
      if (allocated(self%zbotcp)) deallocate(self%zbotcp)
      if (allocated(self%layer))  deallocate(self%layer)
      if (allocated(self%inpola)) deallocate(self%inpola)
      if (allocated(self%inpolb)) deallocate(self%inpolb)
      if (allocated(self%botcom)) deallocate(self%botcom)
      if (allocated(self%nod1lay)) deallocate(self%nod1lay)
      allocate(self%dz(macp))
      allocate(self%z(macp))
      allocate(self%disnod(macp+1))
      allocate(self%ztopcp(macp))
      allocate(self%zbotcp(macp))
      allocate(self%layer(macp))
      allocate(self%inpola(macp));   self%inpola  = 0.0d0
      allocate(self%inpolb(macp));   self%inpolb  = 0.0d0
      allocate(self%botcom(maho));   self%botcom  = 0
      allocate(self%nod1lay(maho));  self%nod1lay = 0
   end subroutine allocate_grid_arrays

   !> Nodal positions and inter-node distances; layer of each node.
   subroutine place_nodes(self, soil_config)
      type(mesh_state_t),  intent(inout) :: self
      type(soil_config_t), intent(in)    :: soil_config

      integer :: i, j, node, nsublay
      real(real64) :: hcomp_i

      nsublay = 0
      if (allocated(soil_config%sublay)) nsublay = size(soil_config%sublay)

      ! Nodal positions and inter-node distances; layer of each node.
      node = 0
      do i = 1, nsublay
         hcomp_i = soil_config%hsublay(i) / dble(soil_config%ncomp(i))
         do j = 1, soil_config%ncomp(i)
            node = node + 1
            self%dz(node) = hcomp_i
            if (node .eq. 1) then
               self%z(node)      = -0.5d0 * self%dz(node)
               self%disnod(node) = -self%z(node)
               self%layer(node)  = soil_config%isoillay(i)
            else
               self%z(node)      = self%z(node-1) - 0.5d0*(self%dz(node-1) + self%dz(node))
               self%disnod(node) = self%z(node-1) - self%z(node)
               self%layer(node)  = soil_config%isoillay(i)
            end if
         end do
      end do
      self%numnod = node
      self%disnod(self%numnod + 1) = 0.5d0 * self%dz(self%numnod)
   end subroutine place_nodes

   !> Top/bottom depths per compartment.
   subroutine compute_compartment_depths(self)
      type(mesh_state_t), intent(inout) :: self

      integer :: i

      ! Top/bottom depths per compartment.
      do i = 1, self%numnod
         if (i == 1) then
            self%ztopcp(i) = 0.0d0
            self%zbotcp(i) = -self%dz(i)
         else
            self%ztopcp(i) = self%zbotcp(i-1)
            self%zbotcp(i) = self%zbotcp(i-1) - self%dz(i)
         end if
      end do
   end subroutine compute_compartment_depths

   !> Bottom compartment and first node of each soil layer.
   subroutine map_layers_to_nodes(self)
      type(mesh_state_t), intent(inout) :: self

      integer :: node, layold, lay

      ! Bottom compartment of each soil layer.
      layold = 1
      do node = 1, self%numnod
         if (self%layer(node) .gt. layold) then
            self%botcom(layold) = node - 1
            layold = layold + 1
         end if
      end do
      self%numlay = layold
      self%botcom(self%numlay) = self%numnod

      ! First node of each soil layer.
      do lay = 1, self%numlay
         node = 1
         do while (self%layer(node) .ne. lay)
            node = node + 1
         end do
         self%nod1lay(lay) = node
      end do
   end subroutine map_layers_to_nodes

   !> Linear-interpolation weights between nodes (mesh-derived).
   subroutine compute_interpolation_weights(self)
      type(mesh_state_t), intent(inout) :: self

      integer :: node

      ! Linear-interpolation weights between nodes (mesh-derived).
      self%inpolb(1) = 0.5d0*self%dz(1)/self%disnod(2)
      do node = 2, self%numnod - 1
         self%inpola(node) = 0.5d0*self%dz(node)/self%disnod(node)
         self%inpolb(node) = 0.5d0*self%dz(node)/self%disnod(node+1)
      end do
      self%inpola(self%numnod) = 0.5d0*self%dz(self%numnod)/self%disnod(self%numnod)
   end subroutine compute_interpolation_weights

end module mesh_state_mod
