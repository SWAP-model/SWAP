!> Surface Water State Module
!!
!! Manages surface-water levels, storage, management controls, and lookup
!! tables used by the extended drainage routines.
!!
!! Key components:
!! - Configuration and management parameters
!! - Current surface-water state variables
!! - Cumulative and intermediate flux tracking
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module surfacewater_state_mod
	implicit none
	private

	public :: surfacewater_state_t
	public :: surfacewater_state_init
	public :: surfacewater_state_finalize
	public :: surfacewater_state_reset_cumulative
	public :: surfacewater_state_reset_intermediate

	!> Surface-water state container.
	type :: surfacewater_state_t
		! Water levels
		real(8) :: wlp = 0.0d0                 !! Primary water level
		real(8) :: wls = 0.0d0                 !! Secondary water level
		real(8) :: wlsold = 0.0d0              !! Previous secondary water level
		real(8) :: wlstar = 0.0d0              !! Target water level
		real(8) :: hwlman = 0.0d0              !! Managed water level
		real(8) :: vtair = 0.0d0               !! Total air volume in profile
		real(8), allocatable :: wlsbak(:)      !! Last 4 water levels for oscillation check

		! Storage
		real(8) :: swst = 0.0d0                !! Surface-water storage
		real(8) :: swstini = 0.0d0             !! Initial surface-water storage

		! Fluxes
		real(8) :: qdrd = 0.0d0                !! Drainage discharge to secondary system
		real(8) :: cqdrd = 0.0d0               !! Cumulative drainage discharge
		real(8) :: cwsupp = 0.0d0              !! Cumulative water supply
		real(8) :: cwout = 0.0d0               !! Cumulative water outflow
		real(8) :: runots = 0.0d0              !! Runoff during current time step
		real(8) :: QRapDra = 0.0d0             !! Rapid drainage flux

		! Management parameters
		integer :: imper = 0                   !! Current management period
		integer :: nmper = 0                   !! Number of management periods
		integer :: numadj = 0                  !! Number of level adjustments
		integer :: swsrf = 0                   !! Surface-water switch
		integer :: swsec = 0                   !! Secondary water switch
		integer :: swqhr = 0                   !! Q-H relation switch
		real(8) :: osswlm = 0.0d0              !! Oscillation tolerance
		real(8), allocatable :: impend(:)      !! Management period end dates
		integer, allocatable :: swman(:)       !! Management type per period
		real(8), allocatable :: hbweir(:)      !! Weir crest level per period
		real(8), allocatable :: wldip(:)       !! Dip below target for supply
		real(8), allocatable :: alphaw(:)      !! Weir discharge coefficient a
		real(8), allocatable :: betaw(:)       !! Weir discharge coefficient b
		real(8), allocatable :: wscap(:)       !! Maximum supply capacity per period
		real(8), allocatable :: dropr(:)       !! Maximum drop rate per period
		real(8), allocatable :: intwl(:)       !! Adjustment interval per period
		integer, allocatable :: nphase(:)      !! Number of phases per period
		real(8), allocatable :: gwlcrit(:,:)   !! Critical GWL per period/phase
		real(8), allocatable :: wlsman(:,:)    !! Target level per period/phase
		real(8), allocatable :: vcrit(:,:)     !! Critical air volume per period/phase
		real(8), allocatable :: hcrit(:,:)     !! Critical pressure head per period/phase
		integer, allocatable :: nodhd(:)       !! Node index for head criterion per period

		! Lookup tables
		real(8), allocatable :: wlstab(:)      !! Secondary water-level table
		real(8), allocatable :: wlptab(:)      !! Primary water-level table
		real(8), allocatable :: sttab(:,:)     !! Storage-level lookup table
		real(8), allocatable :: owltab(:,:)    !! Open-water level tables per drain level
		real(8), allocatable :: qqhtab(:,:)    !! Q-H table per management period

		! Flags
		logical :: flsurfacewater = .false.    !! Surface-water simulation active
		logical :: overfl = .false.            !! Overflow occurred
		logical :: fldecdt = .false.           !! Decrease timestep flag
		logical :: fldtmin = .false.           !! At minimum timestep flag
	end type surfacewater_state_t

contains

	!> Initialize the surface-water state.
	!!
	!! Allocates dynamic arrays and resets all fields to canonical defaults.
	!!
	!! @param[inout] state Surface-water state container
	!! @param[in] nmper Number of management periods
	!! @param[in] mamte Maximum number of phase entries per management period
	!! @param[in] maowl Maximum number of open-water-level table entries
	!! @param[in] nrlevs Number of drainage levels
	subroutine surfacewater_state_init(state, nmper, mamte, maowl, nrlevs)
		type(surfacewater_state_t), intent(inout) :: state
		integer, intent(in) :: nmper
		integer, intent(in) :: mamte
		integer, intent(in) :: maowl
		integer, intent(in) :: nrlevs

		allocate(state%wlsbak(4))

		allocate(state%impend(nmper))
		allocate(state%swman(nmper))
		allocate(state%hbweir(nmper))
		allocate(state%wldip(nmper))
		allocate(state%alphaw(nmper))
		allocate(state%betaw(nmper))
		allocate(state%wscap(nmper))
		allocate(state%dropr(nmper))
		allocate(state%intwl(nmper))
		allocate(state%nphase(nmper))
		allocate(state%nodhd(nmper))
		allocate(state%gwlcrit(nmper, mamte))
		allocate(state%wlsman(nmper, mamte))
		allocate(state%vcrit(nmper, mamte))
		allocate(state%hcrit(nmper, mamte))

		allocate(state%wlstab(2*maowl))
		allocate(state%wlptab(2*maowl))
		allocate(state%sttab(22, 2))
		allocate(state%owltab(nrlevs, 2*maowl))
		allocate(state%qqhtab(nmper, 22))

		call surfacewater_state_reset_all(state)
		state%nmper = nmper

	end subroutine surfacewater_state_init

	!> Finalize the surface-water state.
	!!
	!! Deallocates all dynamic arrays and resets scalar fields to canonical
	!! defaults to avoid carry-over between model instances.
	!!
	!! @param[inout] state Surface-water state container
	subroutine surfacewater_state_finalize(state)
		type(surfacewater_state_t), intent(inout) :: state

		if (allocated(state%wlsbak)) deallocate(state%wlsbak)
		if (allocated(state%impend)) deallocate(state%impend)
		if (allocated(state%swman)) deallocate(state%swman)
		if (allocated(state%hbweir)) deallocate(state%hbweir)
		if (allocated(state%wldip)) deallocate(state%wldip)
		if (allocated(state%alphaw)) deallocate(state%alphaw)
		if (allocated(state%betaw)) deallocate(state%betaw)
		if (allocated(state%wscap)) deallocate(state%wscap)
		if (allocated(state%dropr)) deallocate(state%dropr)
		if (allocated(state%intwl)) deallocate(state%intwl)
		if (allocated(state%nphase)) deallocate(state%nphase)
		if (allocated(state%nodhd)) deallocate(state%nodhd)
		if (allocated(state%gwlcrit)) deallocate(state%gwlcrit)
		if (allocated(state%wlsman)) deallocate(state%wlsman)
		if (allocated(state%vcrit)) deallocate(state%vcrit)
		if (allocated(state%hcrit)) deallocate(state%hcrit)
		if (allocated(state%wlstab)) deallocate(state%wlstab)
		if (allocated(state%wlptab)) deallocate(state%wlptab)
		if (allocated(state%sttab)) deallocate(state%sttab)
		if (allocated(state%owltab)) deallocate(state%owltab)
		if (allocated(state%qqhtab)) deallocate(state%qqhtab)

		call surfacewater_state_reset_all(state)

	end subroutine surfacewater_state_finalize

	!> Reset all surface-water fields to canonical defaults.
	!!
	!! Resets scalar fields and zeroes allocatable arrays if they are already
	!! allocated.
	!!
	!! @param[inout] state Surface-water state container
	subroutine surfacewater_state_reset_all(state)
		type(surfacewater_state_t), intent(inout) :: state

		state%wlp = 0.0d0
		state%wls = 0.0d0
		state%wlsold = 0.0d0
		state%wlstar = 0.0d0
		state%hwlman = 0.0d0
		state%vtair = 0.0d0
		if (allocated(state%wlsbak)) state%wlsbak = 0.0d0

		state%swst = 0.0d0
		state%swstini = 0.0d0
		state%qdrd = 0.0d0
		state%QRapDra = 0.0d0

		state%imper = 0
		state%nmper = 0
		state%numadj = 0
		state%swsrf = 0
		state%swsec = 0
		state%swqhr = 0
		state%osswlm = 0.0d0

		if (allocated(state%impend)) state%impend = 0.0d0
		if (allocated(state%swman)) state%swman = 0
		if (allocated(state%hbweir)) state%hbweir = 0.0d0
		if (allocated(state%wldip)) state%wldip = 0.0d0
		if (allocated(state%alphaw)) state%alphaw = 0.0d0
		if (allocated(state%betaw)) state%betaw = 0.0d0
		if (allocated(state%wscap)) state%wscap = 0.0d0
		if (allocated(state%dropr)) state%dropr = 0.0d0
		if (allocated(state%intwl)) state%intwl = 0.0d0
		if (allocated(state%nphase)) state%nphase = 0
		if (allocated(state%nodhd)) state%nodhd = 0
		if (allocated(state%gwlcrit)) state%gwlcrit = 0.0d0
		if (allocated(state%wlsman)) state%wlsman = 0.0d0
		if (allocated(state%vcrit)) state%vcrit = 0.0d0
		if (allocated(state%hcrit)) state%hcrit = 0.0d0

		if (allocated(state%wlstab)) state%wlstab = 0.0d0
		if (allocated(state%wlptab)) state%wlptab = 0.0d0
		if (allocated(state%sttab)) state%sttab = 0.0d0
		if (allocated(state%owltab)) state%owltab = 0.0d0
		if (allocated(state%qqhtab)) state%qqhtab = 0.0d0

		call surfacewater_state_reset_cumulative(state)
		call surfacewater_state_reset_intermediate(state)

		state%flsurfacewater = .false.
		state%overfl = .false.
		state%fldecdt = .false.
		state%fldtmin = .false.

	end subroutine surfacewater_state_reset_all

	!> Reset cumulative surface-water accumulators.
	!!
	!! @param[inout] state Surface-water state container
	subroutine surfacewater_state_reset_cumulative(state)
		type(surfacewater_state_t), intent(inout) :: state

		state%cqdrd = 0.0d0
		state%cwsupp = 0.0d0
		state%cwout = 0.0d0

	end subroutine surfacewater_state_reset_cumulative

	!> Reset intermediate surface-water accumulators.
	!!
	!! @param[inout] state Surface-water state container
	subroutine surfacewater_state_reset_intermediate(state)
		type(surfacewater_state_t), intent(inout) :: state

		state%runots = 0.0d0

	end subroutine surfacewater_state_reset_intermediate

end module surfacewater_state_mod
