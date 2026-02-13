!> Solute State Module
!!
!! Manages solute concentrations, transport parameters, cumulative and
!! intermediate mass balances, and age-tracer state.
!!
!! Key components:
!! - Solute configuration and transport switches
!! - Node/layer solute and sorption state
!! - Cumulative and intermediate solute mass tracking
!! - Age tracer state variables
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module solute_state_mod
	use swap_array_dimensions, only: MABBC
	implicit none
	private

	public :: solute_state_t
	public :: solute_state_init
	public :: solute_state_finalize
	public :: solute_state_reset_cumulative
	public :: solute_state_reset_intermediate

	!> Solute state container.
	type :: solute_state_t
		! Configuration switches
		integer :: swsolu = 0                  !! Switch for solute transport simulation
		integer :: swsp = 0                    !! Switch for sorption simulation
		integer :: swbr = 0                    !! Switch for breakthrough curve
		integer :: swbotbc = 1                 !! Switch for bottom BC
		integer :: nconc = 0                   !! Number of initial concentrations

		! Concentrations
		real(8), allocatable :: cml(:)         !! Mobile concentration (mg/cm3)
		real(8), allocatable :: cmsy(:)        !! Total concentration
		real(8) :: cpond = 0.0d0               !! Ponding concentration
		real(8) :: csurf = 0.0d0               !! Surface solute amount
		real(8) :: cdrain = 0.0d0              !! Drainage concentration
		real(8) :: cseep = 0.0d0               !! Seepage concentration
		real(8) :: cpre = 0.0d0                !! Precipitation concentration
		real(8) :: cirr = 0.0d0                !! Irrigation concentration
		real(8) :: cref = 1.0d0                !! Reference concentration for Freundlich

		! Cumulative amounts
		real(8) :: sampro = 0.0d0              !! Total in profile
		real(8) :: samini = 0.0d0              !! Initial in profile
		real(8) :: sqbot = 0.0d0               !! Through bottom
		real(8) :: sqdra = 0.0d0               !! To drainage
		real(8) :: sqprec = 0.0d0              !! In precipitation
		real(8) :: sqirrig = 0.0d0             !! In irrigation
		real(8) :: sqsur = 0.0d0               !! To surface water
		real(8) :: sqrap = 0.0d0               !! In rapid drainage
		real(8) :: dectot = 0.0d0              !! Decomposition
		real(8) :: rottot = 0.0d0              !! Root extraction
		real(8) :: solbal = 0.0d0              !! Balance error

		! Intermediate amounts
		real(8) :: imsqbot = 0.0d0
		real(8) :: imsqdra = 0.0d0
		real(8) :: imsqprec = 0.0d0
		real(8) :: imsqirrig = 0.0d0
		real(8) :: imdectot = 0.0d0
		real(8) :: imrottot = 0.0d0
		real(8) :: isqbot = 0.0d0
		real(8) :: isqtop = 0.0d0

		! Macropore solute
		real(8) :: samcra = 0.0d0              !! In cracks

		! Age tracer
		real(8) :: AgeGwl1m = 0.0d0            !! GW age upper 1m sat zone
		real(8) :: icAgeBot = 0.0d0
		real(8) :: icAgeRot = 0.0d0
		real(8) :: icAgeSur = 0.0d0
		real(8), allocatable :: icAgeDra(:)

		! Age tracer boundary/pond state
		real(8) :: Ageirr = 0.0d0              !! Age of irrigation water
		real(8) :: Agedrain = 0.0d0            !! Age of drainage water
		real(8) :: Agepre = 0.0d0              !! Age of precipitation
		real(8) :: Agepond = 0.0d0             !! Age of ponding water
		real(8) :: Agepondm1 = 0.0d0           !! Age of ponding (prev timestep)
		real(8) :: icAgetopupw = 0.0d0         !! Incremental age leaving top
		real(8) :: icAgetopdwn = 0.0d0         !! Incremental age entering top
		real(8) :: ArMpSs = 0.0d0              !! Area fraction macropores at surface

		! Transport parameters
		real(8) :: ddif = 0.0d0                !! Molecular diffusion coefficient
		real(8) :: frexp = 1.0d0               !! Freundlich exponent
		real(8) :: tscf = 1.0d0                !! Relative uptake by roots
		real(8) :: dtsolu = 0.0d0              !! Max time step for solute

		! Decomposition parameters
		real(8) :: gampar = 0.0d0              !! Temp reduction factor
		real(8) :: bexp = 0.7d0                !! Dryness exponent
		real(8) :: rtheta = 0.01d0             !! Min theta for decomposition
		real(8) :: decsat = 0.0d0              !! Decomposition in aquifer

		! Aquifer parameters for breakthrough
		real(8) :: daquif = 0.0d0              !! Aquifer thickness
		real(8) :: poros = 0.0d0               !! Aquifer porosity
		real(8) :: kfsat = 0.0d0               !! Adsorption in aquifer

		! Salt stress parameters
		real(8) :: salthead = 0.0d0            !! Salt to osmotic head
		real(8) :: saltmax = 0.0d0             !! Threshold concentration
		real(8) :: saltslope = 0.0d0           !! Uptake decline

		! Parameters per layer
		real(8), allocatable :: ldis(:)        !! Dispersion length
		real(8), allocatable :: kf(:)          !! Freundlich coefficient
		real(8), allocatable :: decpot(:)      !! Potential decomposition rate
		real(8), allocatable :: fdepth(:)      !! Depth reduction factor

		! Tables
		real(8), allocatable :: cseeptab(:)    !! Seepage concentration table
		real(8), allocatable :: zc(:)          !! Depths for initial concentrations

		! Flags
		logical :: flsolute = .false.
		logical :: flAgeTracer = .false.
	end type solute_state_t

contains

	!> Initialize the solute state.
	!!
	!! Allocates dynamic arrays and resets all fields to canonical defaults.
	!!
	!! @param[inout] state Solute state container
	!! @param[in] numnod Number of soil compartments
	!! @param[in] numlay Number of soil layers
	!! @param[in] nrlevs Number of drainage levels
	!! @param[in] ncseep Optional cseeptab size (default: 2*MABBC)
	subroutine solute_state_init(state, numnod, numlay, nrlevs, ncseep)
		type(solute_state_t), intent(inout) :: state
		integer, intent(in) :: numnod
		integer, intent(in) :: numlay
		integer, intent(in) :: nrlevs
		integer, intent(in), optional :: ncseep

		integer :: ncs

		ncs = 2 * MABBC
		if (present(ncseep)) ncs = ncseep

		allocate(state%cml(numnod))
		allocate(state%cmsy(numnod))
		allocate(state%ldis(numlay))
		allocate(state%kf(numlay))
		allocate(state%decpot(numlay))
		allocate(state%fdepth(numlay))
		allocate(state%icAgeDra(nrlevs))
		allocate(state%cseeptab(ncs))
		allocate(state%zc(numnod))

		call solute_state_reset_all(state)
		state%fdepth = 1.0d0

	end subroutine solute_state_init

	!> Finalize the solute state.
	!!
	!! Deallocates all dynamic arrays and resets scalar fields to canonical
	!! defaults to avoid carry-over between model instances.
	!!
	!! @param[inout] state Solute state container
	subroutine solute_state_finalize(state)
		type(solute_state_t), intent(inout) :: state

		if (allocated(state%cml)) deallocate(state%cml)
		if (allocated(state%cmsy)) deallocate(state%cmsy)
		if (allocated(state%icAgeDra)) deallocate(state%icAgeDra)
		if (allocated(state%ldis)) deallocate(state%ldis)
		if (allocated(state%kf)) deallocate(state%kf)
		if (allocated(state%decpot)) deallocate(state%decpot)
		if (allocated(state%fdepth)) deallocate(state%fdepth)
		if (allocated(state%cseeptab)) deallocate(state%cseeptab)
		if (allocated(state%zc)) deallocate(state%zc)

		call solute_state_reset_all(state)

	end subroutine solute_state_finalize

	!> Reset all solute fields to canonical defaults.
	!!
	!! Resets scalar fields and zeroes allocatable arrays if they are already
	!! allocated.
	!!
	!! @param[inout] state Solute state container
	subroutine solute_state_reset_all(state)
		type(solute_state_t), intent(inout) :: state

		state%swsolu = 0
		state%swsp = 0
		state%swbr = 0
		state%swbotbc = 1
		state%nconc = 0

		if (allocated(state%cml)) state%cml = 0.0d0
		if (allocated(state%cmsy)) state%cmsy = 0.0d0
		state%cpond = 0.0d0
		state%csurf = 0.0d0
		state%cdrain = 0.0d0
		state%cseep = 0.0d0
		state%cpre = 0.0d0
		state%cirr = 0.0d0
		state%cref = 1.0d0

		call solute_state_reset_cumulative(state)
		call solute_state_reset_intermediate(state)

		state%samcra = 0.0d0

		state%AgeGwl1m = 0.0d0
		state%icAgeBot = 0.0d0
		state%icAgeRot = 0.0d0
		state%icAgeSur = 0.0d0
		if (allocated(state%icAgeDra)) state%icAgeDra = 0.0d0

		state%Ageirr = 0.0d0
		state%Agedrain = 0.0d0
		state%Agepre = 0.0d0
		state%Agepond = 0.0d0
		state%Agepondm1 = 0.0d0
		state%icAgetopupw = 0.0d0
		state%icAgetopdwn = 0.0d0
		state%ArMpSs = 0.0d0

		state%ddif = 0.0d0
		state%frexp = 1.0d0
		state%tscf = 1.0d0
		state%dtsolu = 0.0d0

		state%gampar = 0.0d0
		state%bexp = 0.7d0
		state%rtheta = 0.01d0
		state%decsat = 0.0d0

		state%daquif = 0.0d0
		state%poros = 0.0d0
		state%kfsat = 0.0d0

		state%salthead = 0.0d0
		state%saltmax = 0.0d0
		state%saltslope = 0.0d0

		if (allocated(state%ldis)) state%ldis = 0.0d0
		if (allocated(state%kf)) state%kf = 0.0d0
		if (allocated(state%decpot)) state%decpot = 0.0d0
		if (allocated(state%fdepth)) state%fdepth = 0.0d0

		if (allocated(state%cseeptab)) state%cseeptab = 0.0d0
		if (allocated(state%zc)) state%zc = 0.0d0

		state%flsolute = .false.
		state%flAgeTracer = .false.

	end subroutine solute_state_reset_all

	!> Reset cumulative solute mass accumulators.
	!!
	!! @param[inout] state Solute state container
	subroutine solute_state_reset_cumulative(state)
		type(solute_state_t), intent(inout) :: state

		state%sampro = 0.0d0
		state%samini = 0.0d0
		state%sqbot = 0.0d0
		state%sqdra = 0.0d0
		state%sqprec = 0.0d0
		state%sqirrig = 0.0d0
		state%sqsur = 0.0d0
		state%sqrap = 0.0d0
		state%dectot = 0.0d0
		state%rottot = 0.0d0
		state%solbal = 0.0d0

	end subroutine solute_state_reset_cumulative

	!> Reset intermediate solute accumulators.
	!!
	!! @param[inout] state Solute state container
	subroutine solute_state_reset_intermediate(state)
		type(solute_state_t), intent(inout) :: state

		state%imsqbot = 0.0d0
		state%imsqdra = 0.0d0
		state%imsqprec = 0.0d0
		state%imsqirrig = 0.0d0
		state%imdectot = 0.0d0
		state%imrottot = 0.0d0
		state%isqbot = 0.0d0
		state%isqtop = 0.0d0

	end subroutine solute_state_reset_intermediate

end module solute_state_mod
