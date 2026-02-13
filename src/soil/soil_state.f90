!> Soil State Module
!!
!! This module defines the soil water state for SWAP.
!! Manages soil hydraulic state variables, discretization geometry,
!! groundwater and ponding state, and cumulative/intermediate flux tracking.
!!
!! Key components:
!! - Soil hydraulic and storage state per node and layer
!! - Geometry/discretization metadata
!! - Groundwater and boundary-related soil state
!! - Cumulative and intermediate water flux tracking
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module soil_state_mod
	implicit none
	private

	public :: soil_state_t
	public :: soil_state_init
	public :: soil_state_finalize
	public :: soil_state_reset_cumulative
	public :: soil_state_reset_intermediate

	type :: soil_state_t
		! Primary state variables (node-based)
		real(8), allocatable :: h(:)           !! Pressure head (cm)
		real(8), allocatable :: hm1(:)         !! Pressure head at previous time
		real(8), allocatable :: theta(:)       !! Volumetric water content (-)
		real(8), allocatable :: thetm1(:)      !! Water content at previous time
		real(8), allocatable :: k(:)           !! Hydraulic conductivity (cm/d)
		real(8), allocatable :: kmean(:)       !! Mean K at interfaces
		real(8), allocatable :: dimoca(:)      !! Differential moisture capacity

		! Geometry
		real(8), allocatable :: z(:)           !! Node depths (cm)
		real(8), allocatable :: dz(:)          !! Compartment thickness (cm)
		real(8), allocatable :: ztopcp(:)      !! Top of compartment
		real(8), allocatable :: zbotcp(:)      !! Bottom of compartment
		real(8), allocatable :: disnod(:)      !! Distance between nodes

		! Soil properties per node
		real(8), allocatable :: thetar(:)      !! Residual water content
		real(8), allocatable :: thetas(:)      !! Saturated water content
		real(8), allocatable :: hroot(:)       !! Pressure head at root-soil interface
		real(8), allocatable :: twilt(:)       !! Wilting point pressure head
		real(8), allocatable :: rfcp(:)        !! Frost reduction factor

		! Soil layer properties
		real(8), allocatable :: bdens(:)       !! Bulk density (g/cm3)
		real(8), allocatable :: cofani(:)      !! Anisotropy coefficient
		real(8), allocatable :: ksatfit(:)     !! Saturated K fitted
		real(8), allocatable :: ksatexm(:)     !! Saturated K measured
		real(8), allocatable :: thetsl(:)      !! Saturated water content per layer

		! Van Genuchten parameters
		real(8), allocatable :: paramvg(:,:)   !! VG parameters (21,maho)
		real(8), allocatable :: cofgen(:,:)    !! Adjusted VG parameters (21,macp)

		! Fluxes
		real(8), allocatable :: q(:)           !! Water flux between compartments
		real(8), allocatable :: evp(:)         !! Evaporation flux
		real(8), allocatable :: qrot(:)        !! Root extraction flux
		real(8), allocatable :: qpotrot(:)     !! Potential root extraction
		real(8), allocatable :: qssdi(:)       !! Subsurface drip irrigation

		! Stress reduction fluxes
		real(8), allocatable :: qredwet(:)     !! Reduction due to wet conditions
		real(8), allocatable :: qreddry(:)     !! Reduction due to dry conditions
		real(8), allocatable :: qredsol(:)     !! Reduction due to salinity
		real(8), allocatable :: qredfrs(:)     !! Reduction due to frost

		! Intermediate storage
		real(8), allocatable :: inq(:)         !! Intermediate water flux
		real(8), allocatable :: inqrot(:)      !! Intermediate root extraction
		real(8), allocatable :: inqssdi(:)     !! Intermediate SSDI
		real(8), allocatable :: ithetabeg(:)   !! Water content at period start

		! Groundwater
		real(8) :: gwl = 0.0d0                 !! Groundwater level (cm)
		real(8) :: gwlm1 = 0.0d0               !! GWL at previous time
		real(8) :: gwli = 0.0d0                !! Initial GWL
		real(8) :: gwlinp = 0.0d0              !! Prescribed GWL
		real(8) :: pegwl = 0.0d0               !! Perched GWL
		real(8) :: deepgw = 0.0d0              !! Hydraulic head in aquifer
		integer :: nodgwl = 0                  !! Node above GWL
		integer :: npegwl = 0                  !! Node above perched GWL
		integer :: bpegwl = 0                  !! Node at bottom perched GWL

		! Surface water/ponding
		real(8) :: pond = 0.0d0                !! Ponding depth (cm)
		real(8) :: pondm1 = 0.0d0              !! Ponding at previous time
		real(8) :: pondmx = 0.0d0              !! Maximum ponding
		real(8) :: pondini = 0.0d0             !! Initial ponding
		real(8) :: hsurf = 0.0d0               !! Pressure head at surface

		! Bottom boundary
		real(8) :: qbot = 0.0d0                !! Bottom flux
		real(8) :: hbot = 0.0d0                !! Bottom pressure head

		! Top boundary
		real(8) :: qtop = 0.0d0                !! Top flux
		logical :: ftoph = .false.             !! Prescribed pressure head at top

		! Cumulative fluxes
		real(8) :: cqbot = 0.0d0               !! Cumulative bottom flux
		real(8) :: cqbotdo = 0.0d0             !! Cumulative downward bottom
		real(8) :: cqbotup = 0.0d0             !! Cumulative upward bottom
		real(8) :: cqtdo = 0.0d0               !! Cumulative downward top
		real(8) :: cqtup = 0.0d0               !! Cumulative upward top
		real(8) :: cqrot = 0.0d0               !! Cumulative root extraction
		real(8) :: cqdra = 0.0d0               !! Cumulative drainage
		real(8) :: crunoff = 0.0d0             !! Cumulative runoff
		real(8) :: crunon = 0.0d0              !! Cumulative runon

		! Intermediate fluxes
		real(8) :: iqbot = 0.0d0
		real(8) :: iqrot = 0.0d0
		real(8) :: iqdra = 0.0d0
		real(8) :: iruno = 0.0d0
		real(8) :: irunon = 0.0d0
		real(8) :: iqssdi = 0.0d0
		real(8) :: ipondbeg = 0.0d0

		! Storage
		real(8) :: volact = 0.0d0              !! Current water storage
		real(8) :: volini = 0.0d0              !! Initial water storage
		real(8) :: volm1 = 0.0d0               !! Previous water storage
		real(8) :: wbalance = 0.0d0            !! Water balance error

		! Iteration control
		integer :: numbit = 0                  !! Iteration number
		integer :: msteps = 0                  !! Max iteration steps per day
		real(8) :: CritDevh1Cp = 0.0d0         !! Convergence criterion h relative
		real(8) :: CritDevh2Cp = 0.0d0         !! Convergence criterion h absolute
		real(8) :: CritDevMasBal = 0.0d0       !! Max water balance error
		real(8) :: gwlconv = 0.0d0             !! GWL convergence criterion

		! Discretization info
		integer :: numnod = 0                  !! Number of nodes
		integer :: numlay = 0                  !! Number of soil layers
		integer :: nsublay = 0                 !! Number of sublayers
		integer, allocatable :: layer(:)       !! Layer number per compartment
		integer, allocatable :: botcom(:)      !! Bottom compartment per layer
		integer, allocatable :: nod1lay(:)     !! First node per layer

		! Hysteresis
		integer :: swhyst = 0                  !! Hysteresis switch
		integer, allocatable :: indeks(:)      !! Wetting/drying curve index
		real(8) :: tau = 0.0d0                 !! Min h difference for transition

		! Soil evaporation reduction
		real(8) :: saev = 0.0d0                !! Cumulative actual E (Boesten)
		real(8) :: spev = 0.0d0                !! Cumulative potential E (Boesten)
		real(8) :: ldwet = 0.0d0               !! Dry period length (Black)
		real(8) :: cofred = 0.0d0              !! Reduction coefficient

		! Headcalc iteration tracking (from headcalc.f90 SAVE variables)
		logical :: flwarn_hc = .true.          !! Warning flag for headcalc
		integer :: iwarn_hc = 0                !! Warning counter
		integer :: nstep_hc = 0                !! Step counter in headcalc

		! Flags
		logical :: FlRunoff = .false.
		logical :: fldrain = .false.
		logical :: flrunon = .false.
		logical :: fllowgwl = .false.
	end type soil_state_t

contains

	!> Initialize the soil state.
	!!
	!! Allocates node-based and layer-based arrays and sets all values to
	!! canonical defaults.
	!!
	!! @param[inout] state Soil state container
	!! @param[in] numnod Number of soil nodes
	!! @param[in] numlay Number of soil layers
	subroutine soil_state_init(state, numnod, numlay)
		type(soil_state_t), intent(inout) :: state
		integer, intent(in) :: numnod
		integer, intent(in) :: numlay

		allocate(state%h(numnod))
		allocate(state%hm1(numnod))
		allocate(state%theta(numnod))
		allocate(state%thetm1(numnod))
		allocate(state%k(numnod + 1))
		allocate(state%kmean(numnod + 1))
		allocate(state%dimoca(numnod))
		allocate(state%z(numnod))
		allocate(state%dz(numnod))
		allocate(state%ztopcp(numnod))
		allocate(state%zbotcp(numnod))
		allocate(state%disnod(numnod + 1))
		allocate(state%thetar(numnod))
		allocate(state%thetas(numnod))
		allocate(state%hroot(numnod))
		allocate(state%twilt(numnod))
		allocate(state%rfcp(numnod))
		allocate(state%q(numnod + 1))
		allocate(state%evp(numnod))
		allocate(state%qrot(numnod))
		allocate(state%qpotrot(numnod))
		allocate(state%qssdi(numnod))
		allocate(state%qredwet(numnod))
		allocate(state%qreddry(numnod))
		allocate(state%qredsol(numnod))
		allocate(state%qredfrs(numnod))
		allocate(state%inq(numnod + 1))
		allocate(state%inqrot(numnod))
		allocate(state%inqssdi(numnod))
		allocate(state%ithetabeg(numnod))
		allocate(state%layer(numnod))
		allocate(state%indeks(numnod))

		allocate(state%bdens(numlay))
		allocate(state%cofani(numlay))
		allocate(state%ksatfit(numlay))
		allocate(state%ksatexm(numlay))
		allocate(state%thetsl(numlay))
		allocate(state%botcom(numlay))
		allocate(state%nod1lay(numlay))
		allocate(state%paramvg(21, numlay))
		allocate(state%cofgen(21, numnod))

		call soil_state_reset_all(state)
		state%numnod = numnod
		state%numlay = numlay

	end subroutine soil_state_init

	!> Finalize the soil state.
	!!
	!! Deallocates all dynamic arrays and resets scalar values to canonical
	!! defaults to prevent state carry-over between model instances.
	!!
	!! @param[inout] state Soil state container
	subroutine soil_state_finalize(state)
		type(soil_state_t), intent(inout) :: state

		if (allocated(state%h)) deallocate(state%h)
		if (allocated(state%hm1)) deallocate(state%hm1)
		if (allocated(state%theta)) deallocate(state%theta)
		if (allocated(state%thetm1)) deallocate(state%thetm1)
		if (allocated(state%k)) deallocate(state%k)
		if (allocated(state%kmean)) deallocate(state%kmean)
		if (allocated(state%dimoca)) deallocate(state%dimoca)
		if (allocated(state%z)) deallocate(state%z)
		if (allocated(state%dz)) deallocate(state%dz)
		if (allocated(state%ztopcp)) deallocate(state%ztopcp)
		if (allocated(state%zbotcp)) deallocate(state%zbotcp)
		if (allocated(state%disnod)) deallocate(state%disnod)
		if (allocated(state%thetar)) deallocate(state%thetar)
		if (allocated(state%thetas)) deallocate(state%thetas)
		if (allocated(state%hroot)) deallocate(state%hroot)
		if (allocated(state%twilt)) deallocate(state%twilt)
		if (allocated(state%rfcp)) deallocate(state%rfcp)
		if (allocated(state%bdens)) deallocate(state%bdens)
		if (allocated(state%cofani)) deallocate(state%cofani)
		if (allocated(state%ksatfit)) deallocate(state%ksatfit)
		if (allocated(state%ksatexm)) deallocate(state%ksatexm)
		if (allocated(state%thetsl)) deallocate(state%thetsl)
		if (allocated(state%paramvg)) deallocate(state%paramvg)
		if (allocated(state%cofgen)) deallocate(state%cofgen)
		if (allocated(state%q)) deallocate(state%q)
		if (allocated(state%evp)) deallocate(state%evp)
		if (allocated(state%qrot)) deallocate(state%qrot)
		if (allocated(state%qpotrot)) deallocate(state%qpotrot)
		if (allocated(state%qssdi)) deallocate(state%qssdi)
		if (allocated(state%qredwet)) deallocate(state%qredwet)
		if (allocated(state%qreddry)) deallocate(state%qreddry)
		if (allocated(state%qredsol)) deallocate(state%qredsol)
		if (allocated(state%qredfrs)) deallocate(state%qredfrs)
		if (allocated(state%inq)) deallocate(state%inq)
		if (allocated(state%inqrot)) deallocate(state%inqrot)
		if (allocated(state%inqssdi)) deallocate(state%inqssdi)
		if (allocated(state%ithetabeg)) deallocate(state%ithetabeg)
		if (allocated(state%layer)) deallocate(state%layer)
		if (allocated(state%botcom)) deallocate(state%botcom)
		if (allocated(state%nod1lay)) deallocate(state%nod1lay)
		if (allocated(state%indeks)) deallocate(state%indeks)

		call soil_state_reset_all(state)

	end subroutine soil_state_finalize

	!> Reset all soil state fields to canonical defaults.
	!!
	!! Used by both initialization and finalization. Allocatable arrays are
	!! reset only when currently allocated.
	!!
	!! @param[inout] state Soil state container
	subroutine soil_state_reset_all(state)
		type(soil_state_t), intent(inout) :: state

		if (allocated(state%h)) state%h = 0.0d0
		if (allocated(state%hm1)) state%hm1 = 0.0d0
		if (allocated(state%theta)) state%theta = 0.0d0
		if (allocated(state%thetm1)) state%thetm1 = 0.0d0
		if (allocated(state%k)) state%k = 0.0d0
		if (allocated(state%kmean)) state%kmean = 0.0d0
		if (allocated(state%dimoca)) state%dimoca = 0.0d0
		if (allocated(state%z)) state%z = 0.0d0
		if (allocated(state%dz)) state%dz = 0.0d0
		if (allocated(state%ztopcp)) state%ztopcp = 0.0d0
		if (allocated(state%zbotcp)) state%zbotcp = 0.0d0
		if (allocated(state%disnod)) state%disnod = 0.0d0
		if (allocated(state%thetar)) state%thetar = 0.0d0
		if (allocated(state%thetas)) state%thetas = 0.0d0
		if (allocated(state%hroot)) state%hroot = 0.0d0
		if (allocated(state%twilt)) state%twilt = 0.0d0
		if (allocated(state%rfcp)) state%rfcp = 0.0d0
		if (allocated(state%bdens)) state%bdens = 0.0d0
		if (allocated(state%cofani)) state%cofani = 0.0d0
		if (allocated(state%ksatfit)) state%ksatfit = 0.0d0
		if (allocated(state%ksatexm)) state%ksatexm = 0.0d0
		if (allocated(state%thetsl)) state%thetsl = 0.0d0
		if (allocated(state%paramvg)) state%paramvg = 0.0d0
		if (allocated(state%cofgen)) state%cofgen = 0.0d0
		if (allocated(state%q)) state%q = 0.0d0
		if (allocated(state%evp)) state%evp = 0.0d0
		if (allocated(state%qrot)) state%qrot = 0.0d0
		if (allocated(state%qpotrot)) state%qpotrot = 0.0d0
		if (allocated(state%qssdi)) state%qssdi = 0.0d0
		if (allocated(state%qredwet)) state%qredwet = 0.0d0
		if (allocated(state%qreddry)) state%qreddry = 0.0d0
		if (allocated(state%qredsol)) state%qredsol = 0.0d0
		if (allocated(state%qredfrs)) state%qredfrs = 0.0d0
		if (allocated(state%inq)) state%inq = 0.0d0
		if (allocated(state%inqrot)) state%inqrot = 0.0d0
		if (allocated(state%inqssdi)) state%inqssdi = 0.0d0
		if (allocated(state%ithetabeg)) state%ithetabeg = 0.0d0
		if (allocated(state%layer)) state%layer = 0
		if (allocated(state%botcom)) state%botcom = 0
		if (allocated(state%nod1lay)) state%nod1lay = 0
		if (allocated(state%indeks)) state%indeks = 0

		state%gwl = 0.0d0
		state%gwlm1 = 0.0d0
		state%gwli = 0.0d0
		state%gwlinp = 0.0d0
		state%pegwl = 0.0d0
		state%deepgw = 0.0d0
		state%nodgwl = 0
		state%npegwl = 0
		state%bpegwl = 0

		state%pond = 0.0d0
		state%pondm1 = 0.0d0
		state%pondmx = 0.0d0
		state%pondini = 0.0d0
		state%hsurf = 0.0d0
		state%qbot = 0.0d0
		state%hbot = 0.0d0
		state%qtop = 0.0d0
		state%ftoph = .false.

		call soil_state_reset_cumulative(state)
		call soil_state_reset_intermediate(state)

		state%volact = 0.0d0
		state%volini = 0.0d0
		state%volm1 = 0.0d0
		state%wbalance = 0.0d0

		state%numbit = 0
		state%msteps = 0
		state%CritDevh1Cp = 0.0d0
		state%CritDevh2Cp = 0.0d0
		state%CritDevMasBal = 0.0d0
		state%gwlconv = 0.0d0

		state%numnod = 0
		state%numlay = 0
		state%nsublay = 0

		state%swhyst = 0
		state%tau = 0.0d0

		state%saev = 0.0d0
		state%spev = 0.0d0
		state%ldwet = 0.0d0
		state%cofred = 0.0d0

		state%flwarn_hc = .true.
		state%iwarn_hc = 0
		state%nstep_hc = 0

		state%FlRunoff = .false.
		state%fldrain = .false.
		state%flrunon = .false.
		state%fllowgwl = .false.

	end subroutine soil_state_reset_all

	!> Reset cumulative soil water flux accumulators.
	!!
	!! @param[inout] state Soil state container
	subroutine soil_state_reset_cumulative(state)
		type(soil_state_t), intent(inout) :: state

		state%cqbot = 0.0d0
		state%cqbotdo = 0.0d0
		state%cqbotup = 0.0d0
		state%cqtdo = 0.0d0
		state%cqtup = 0.0d0
		state%cqrot = 0.0d0
		state%cqdra = 0.0d0
		state%crunoff = 0.0d0
		state%crunon = 0.0d0

	end subroutine soil_state_reset_cumulative

	!> Reset intermediate soil water flux accumulators.
	!!
	!! @param[inout] state Soil state container
	subroutine soil_state_reset_intermediate(state)
		type(soil_state_t), intent(inout) :: state

		state%iqbot = 0.0d0
		state%iqrot = 0.0d0
		state%iqdra = 0.0d0
		state%iruno = 0.0d0
		state%irunon = 0.0d0
		state%iqssdi = 0.0d0
		state%ipondbeg = 0.0d0

	end subroutine soil_state_reset_intermediate

end module soil_state_mod
