!> Module for grid and discretization operations
!!
!! This module provides subroutines for calculating grid parameters and
!! converting vertical discretization in soil profiles.
module soilgrid_mod
   use error_mod, only: fatalerr_collected
    implicit none
    private
    public :: calcgrid
contains

      !> Calculate grid parameters for soil compartments
      !!
      !! Computes grid parameters including nodal positions, distances between nodes,
      !! layer assignments, and interpolation coefficients for the soil profile.
      !!
      !! @note
      !! Date: Aug 2004
      !! Purpose: calculate grid parameters
      !! @endnote
      !!
      !! [GR-BH Task 35] signature gained state arg; mesh globals written directly
      !! to state%mesh%X. Legacy globals numlay/botcom/inpola/inpolb/nod1lay kept.
      subroutine calcgrid(state, config)
      use swap_state_mod,        only: swap_state_t
      use swap_config_mod,       only: swap_config_t
      use swap_array_dimensions, only: macp, maho
      implicit none

      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      integer :: i, j, lay, node, layold, nsublay
      real(8) :: hcomp_i
      character(len=200) :: messag
      character(len=11)  :: tmp

      associate (mesh => state%mesh, soil_cfg => config%soil)

         ! sublay-derived counts come from the config arrays.
         nsublay = 0
         if (allocated(soil_cfg%sublay)) nsublay = size(soil_cfg%sublay)

         ! Check correct input of number and height of soil compartments.
         do i = 1, nsublay
            hcomp_i = soil_cfg%hsublay(i) / dble(soil_cfg%ncomp(i))
            if (abs(soil_cfg%ncomp(i)*hcomp_i - soil_cfg%hsublay(i)) .gt. 1.d-5) then
               write(tmp, '(i11)') i
               tmp = adjustl(tmp)
               messag = 'At the soil water section, part 4, at layer '//   &
            &           trim(tmp)//' the height of this soil layer hsublay corresponds'// &
            &           ' not to the product of height and number of compartments'
               call fatalerr_collected('calcgrid', messag)
            end if
         end do

         ! Allocate mesh arrays sized for macp (upper bound; trimmed below).
         if (allocated(mesh%dz))     deallocate(mesh%dz)
         if (allocated(mesh%z))      deallocate(mesh%z)
         if (allocated(mesh%disnod)) deallocate(mesh%disnod)
         if (allocated(mesh%ztopcp)) deallocate(mesh%ztopcp)
         if (allocated(mesh%zbotcp)) deallocate(mesh%zbotcp)
         if (allocated(mesh%layer))  deallocate(mesh%layer)
         if (allocated(mesh%inpola)) deallocate(mesh%inpola)
         if (allocated(mesh%inpolb)) deallocate(mesh%inpolb)
         if (allocated(mesh%botcom)) deallocate(mesh%botcom)
         if (allocated(mesh%nod1lay)) deallocate(mesh%nod1lay)
         allocate(mesh%dz(macp))
         allocate(mesh%z(macp))
         allocate(mesh%disnod(macp+1))
         allocate(mesh%ztopcp(macp))
         allocate(mesh%zbotcp(macp))
         allocate(mesh%layer(macp))
         allocate(mesh%inpola(macp));   mesh%inpola  = 0.0d0
         allocate(mesh%inpolb(macp));   mesh%inpolb  = 0.0d0
         allocate(mesh%botcom(maho));   mesh%botcom  = 0
         allocate(mesh%nod1lay(maho));  mesh%nod1lay = 0

         ! Nodal positions and inter-node distances; layer of each node.
         node = 0
         do i = 1, nsublay
            hcomp_i = soil_cfg%hsublay(i) / dble(soil_cfg%ncomp(i))
            do j = 1, soil_cfg%ncomp(i)
               node = node + 1
               mesh%dz(node) = hcomp_i
               if (node .eq. 1) then
                  mesh%z(node)      = -0.5d0 * mesh%dz(node)
                  mesh%disnod(node) = -mesh%z(node)
                  mesh%layer(node)  = soil_cfg%isoillay(i)
               else
                  mesh%z(node)      = mesh%z(node-1) - 0.5d0*(mesh%dz(node-1) + mesh%dz(node))
                  mesh%disnod(node) = mesh%z(node-1) - mesh%z(node)
                  mesh%layer(node)  = soil_cfg%isoillay(i)
               end if
            end do
         end do
         mesh%numnod = node
         mesh%disnod(mesh%numnod + 1) = 0.5d0 * mesh%dz(mesh%numnod)

         ! Top/bottom depths per compartment.
         do i = 1, mesh%numnod
            if (i == 1) then
               mesh%ztopcp(i) = 0.0d0
               mesh%zbotcp(i) = -mesh%dz(i)
            else
               mesh%ztopcp(i) = mesh%zbotcp(i-1)
               mesh%zbotcp(i) = mesh%zbotcp(i-1) - mesh%dz(i)
            end if
         end do

         ! Bottom compartment of each soil layer.
         layold = 1
         do node = 1, mesh%numnod
            if (mesh%layer(node) .gt. layold) then
               mesh%botcom(layold) = node - 1
               layold = layold + 1
            end if
         end do
         mesh%numlay = layold
         mesh%botcom(mesh%numlay) = mesh%numnod

         ! Linear-interpolation weights between nodes (mesh-derived).
         mesh%inpolb(1) = 0.5d0*mesh%dz(1)/mesh%disnod(2)
         do node = 2, mesh%numnod - 1
            mesh%inpola(node) = 0.5d0*mesh%dz(node)/mesh%disnod(node)
            mesh%inpolb(node) = 0.5d0*mesh%dz(node)/mesh%disnod(node+1)
         end do
         mesh%inpola(mesh%numnod) = 0.5d0*mesh%dz(mesh%numnod)/mesh%disnod(mesh%numnod)

         ! First node of each soil layer.
         do lay = 1, mesh%numlay
            Node = 1
            do while (mesh%layer(Node) .ne. lay)
               Node = Node + 1
            end do
            mesh%nod1lay(lay) = node
         end do

      end associate

      return
      end subroutine calcgrid

    

end module soilgrid_mod