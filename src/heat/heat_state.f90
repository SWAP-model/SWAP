!> Heat State Module
!!
!! This module defines the explicit heat-flow state container used by SWAP.
!! It centralizes allocation, initialization, and finalization of heat
!! simulation data.
!!
!! Key components:
!! - Heat state configuration, temperature profiles, and frost tracking
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module heat_state_mod
	implicit none
	private

	include 'arrays.fi'

	public :: heat_state_t
	public :: heat_state_init
	public :: heat_state_finalize

	!> Heat flow state container.
	!!
	!! Stores switches, boundary parameters, per-compartment thermal
	!! properties, and frost diagnostics.
	type :: heat_state_t
		! Configuration switches
		integer :: swhea = 0                   ! Switch for heat flow simulation
		integer :: swcalt = 1                  ! Method: 1=analytical, 2=numerical
		integer :: swtopbhea = 1               ! Top BC: 1=air temp, 2=measured
		integer :: swbotbhea = 1               ! Bottom BC: 1=zero flux, 2=prescribed
		integer :: swfrost = 0                 ! Switch for frost reduction
		integer :: nheat = 0                   ! Number of initial temperatures

		! Soil temperatures
		real(8), allocatable :: tsoil(:)       ! Temperature per compartment (°C)
		real(8) :: tetop = 0.0d0               ! Top temperature
		real(8) :: tebot = 0.0d0               ! Bottom temperature

		! Thermal properties per compartment
		real(8), allocatable :: heacap(:)      ! Heat capacity (J/cm3/K)
		real(8), allocatable :: heacon(:)      ! Heat conductivity (J/cm/K/d)

		! Frost reduction factor per compartment
		real(8), allocatable :: rfcp(:)        ! Reduction factor for frozen conditions

		! Soil composition per compartment
		real(8), allocatable :: fclay(:)       ! Clay content
		real(8), allocatable :: forg(:)        ! Organic matter content
		real(8), allocatable :: fquartz(:)     ! Sand+silt content

		! Soil composition per layer
		real(8), allocatable :: pclay(:)       ! Clay content per layer
		real(8), allocatable :: psand(:)       ! Sand content per layer
		real(8), allocatable :: psilt(:)       ! Silt content per layer
		real(8), allocatable :: orgmat(:)      ! Organic matter per layer

		! Boundary conditions
		real(8) :: tmean = 0.0d0               ! Mean annual surface temperature
		real(8) :: tampli = 0.0d0              ! Surface temperature amplitude
		real(8) :: timref = 0.0d0              ! Time of max temperature
		real(8) :: ddamp = 0.0d0               ! Damping depth

		! Boundary condition tables
		real(8), allocatable :: tembtab(:)     ! Bottom temperature table
		real(8), allocatable :: temtoptab(:)   ! Top temperature table
		real(8), allocatable :: zh(:)          ! Depths for initial temperatures

		! Frost
		real(8) :: zfrosttop = 0.0d0           ! Top of frost layer
		real(8) :: zfrostbot = 0.0d0           ! Bottom of frost layer
		real(8) :: tfroststa = 0.0d0           ! Frost start temperature
		real(8) :: tfrostend = 0.0d0           ! Frost end temperature
		integer :: nodfrostbot = 0             ! Deepest frost node

		! Flags
		logical :: fltemperature = .false.
	end type heat_state_t

contains

	!> Initialize heat state arrays and defaults.
	!!
	!! Allocates all dynamic arrays and sets canonical values used by
	!! legacy SWAP initialization.
	!!
	!! @param[inout] heat Heat state container
	!! @param[in] numnod Number of soil nodes
	!! @param[in] numlay Number of soil layers
	!!
	!! @note Preserves legacy defaults: `tsoil=10` and `rfcp=1`.
	subroutine heat_state_init(heat, numnod, numlay)
		type(heat_state_t), intent(inout) :: heat
		integer, intent(in) :: numnod
		integer, intent(in) :: numlay

		allocate(heat%tsoil(numnod))
		allocate(heat%heacap(numnod))
		allocate(heat%heacon(numnod))
		allocate(heat%rfcp(numnod))
		allocate(heat%fclay(numnod))
		allocate(heat%forg(numnod))
		allocate(heat%fquartz(numnod))
		allocate(heat%pclay(numlay))
		allocate(heat%psand(numlay))
		allocate(heat%psilt(numlay))
		allocate(heat%orgmat(numlay))
		allocate(heat%tembtab(2 * MABBC))
		allocate(heat%temtoptab(2 * MABBC))
		allocate(heat%zh(numnod))

		call heat_state_reset_all(heat)

		heat%tsoil = 10.0d0
		heat%rfcp = 1.0d0

	end subroutine heat_state_init

	!> Finalize heat state.
	!!
	!! Deallocates dynamic arrays and resets scalar fields.
	!!
	!! @param[inout] heat Heat state container
	subroutine heat_state_finalize(heat)
		type(heat_state_t), intent(inout) :: heat

		if (allocated(heat%tsoil)) deallocate(heat%tsoil)
		if (allocated(heat%heacap)) deallocate(heat%heacap)
		if (allocated(heat%heacon)) deallocate(heat%heacon)
		if (allocated(heat%rfcp)) deallocate(heat%rfcp)
		if (allocated(heat%fclay)) deallocate(heat%fclay)
		if (allocated(heat%forg)) deallocate(heat%forg)
		if (allocated(heat%fquartz)) deallocate(heat%fquartz)
		if (allocated(heat%pclay)) deallocate(heat%pclay)
		if (allocated(heat%psand)) deallocate(heat%psand)
		if (allocated(heat%psilt)) deallocate(heat%psilt)
		if (allocated(heat%orgmat)) deallocate(heat%orgmat)
		if (allocated(heat%tembtab)) deallocate(heat%tembtab)
		if (allocated(heat%temtoptab)) deallocate(heat%temtoptab)
		if (allocated(heat%zh)) deallocate(heat%zh)

		call heat_state_reset_all(heat)

	end subroutine heat_state_finalize


	!> Reset all heat state fields to canonical defaults.
	!!
	!! @param[inout] heat Heat state container
	subroutine heat_state_reset_all(heat)
		type(heat_state_t), intent(inout) :: heat

		heat%swhea = 0
		heat%swcalt = 1
		heat%swtopbhea = 1
		heat%swbotbhea = 1
		heat%swfrost = 0
		heat%nheat = 0
		heat%tetop = 0.0d0
		heat%tebot = 0.0d0
		heat%tmean = 0.0d0
		heat%tampli = 0.0d0
		heat%timref = 0.0d0
		heat%ddamp = 0.0d0
		heat%zfrosttop = 0.0d0
		heat%zfrostbot = 0.0d0
		heat%tfroststa = 0.0d0
		heat%tfrostend = 0.0d0
		heat%nodfrostbot = 0
		heat%fltemperature = .false.

		if (allocated(heat%tsoil)) heat%tsoil = 0.0d0
		if (allocated(heat%heacap)) heat%heacap = 0.0d0
		if (allocated(heat%heacon)) heat%heacon = 0.0d0
		if (allocated(heat%rfcp)) heat%rfcp = 0.0d0
		if (allocated(heat%fclay)) heat%fclay = 0.0d0
		if (allocated(heat%forg)) heat%forg = 0.0d0
		if (allocated(heat%fquartz)) heat%fquartz = 0.0d0
		if (allocated(heat%pclay)) heat%pclay = 0.0d0
		if (allocated(heat%psand)) heat%psand = 0.0d0
		if (allocated(heat%psilt)) heat%psilt = 0.0d0
		if (allocated(heat%orgmat)) heat%orgmat = 0.0d0
		if (allocated(heat%tembtab)) heat%tembtab = 0.0d0
		if (allocated(heat%temtoptab)) heat%temtoptab = 0.0d0
		if (allocated(heat%zh)) heat%zh = 0.0d0

	end subroutine heat_state_reset_all

end module heat_state_mod
