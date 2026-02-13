!> Macropore State Module
!!
!! Defines the explicit macropore state container and lifecycle helpers used by
!! SWAP. This module centralizes macropore configuration, flux bookkeeping,
!! and persistent work arrays that previously relied on hidden `SAVE` state.
!!
!! Key components:
!! - Domain and compartment water storage
!! - Instantaneous, cumulative, and intermediate macropore fluxes
!! - Domain/compartment work arrays formerly stored as local `SAVE`
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module macropore_state_mod
	implicit none
	private

	public :: macropore_state_t
	public :: macropore_state_init
	public :: macropore_state_finalize
	public :: macropore_state_reset_cumulative
	public :: macropore_state_reset_intermediate

	!> Explicit state for macropore flow and storage.
	type :: macropore_state_t
		! Domain water storage
		real(8) :: VlMp = 0.0d0                !! Total macropore volume
		real(8) :: VlMpDm1 = 0.0d0             !! MB domain volume
		real(8) :: VlMpDm2 = 0.0d0             !! IC domain volume
		real(8) :: WaSrDm1 = 0.0d0             !! MB water storage
		real(8) :: WaSrDm2 = 0.0d0             !! IC water storage
		real(8) :: WaSrDm1Ini = 0.0d0          !! Initial MB storage
		real(8) :: WaSrDm2Ini = 0.0d0          !! Initial IC storage
		real(8) :: WaLevDm1 = 0.0d0            !! Water level in MB

		! Fluxes
		real(8) :: QMaPo = 0.0d0               !! Matrix-macropore exchange
		real(8) :: QRapDra = 0.0d0             !! Rapid drainage flux
		real(8) :: QMpLatSs = 0.0d0            !! Lateral inflow at surface
		real(8) :: QInTopLatDm1 = 0.0d0        !! Top lateral inflow MB
		real(8) :: QInTopLatDm2 = 0.0d0        !! Top lateral inflow IC
		real(8) :: QInTopVrtDm1 = 0.0d0        !! Top vertical inflow MB
		real(8) :: QInTopVrtDm2 = 0.0d0        !! Top vertical inflow IC

		! Cumulative fluxes
		real(8) :: cQMpLatSs = 0.0d0
		real(8) :: cQMpOutDrRap = 0.0d0
		real(8) :: cQMpInMtxSatDm1 = 0.0d0
		real(8) :: cQMpInMtxSatDm2 = 0.0d0
		real(8) :: cQMpOutMtxUnsDm1 = 0.0d0
		real(8) :: cQMpOutMtxUnsDm2 = 0.0d0
		real(8) :: cQMpInIntSatDm1 = 0.0d0
		real(8) :: cQMpInIntSatDm2 = 0.0d0
		real(8) :: cQMpInTopLatDm1 = 0.0d0
		real(8) :: cQMpInTopLatDm2 = 0.0d0
		real(8) :: cQMpInTopVrtDm1 = 0.0d0
		real(8) :: cQMpInTopVrtDm2 = 0.0d0
		real(8) :: cQMpOutMtxSatDm1 = 0.0d0
		real(8) :: cQMpOutMtxSatDm2 = 0.0d0

		! Incremental fluxes
		real(8) :: iQMpOutDrRap = 0.0d0
		real(8) :: iQInTopLatDm1 = 0.0d0
		real(8) :: iQInTopLatDm2 = 0.0d0
		real(8) :: iQInTopVrtDm1 = 0.0d0
		real(8) :: iQInTopVrtDm2 = 0.0d0
		real(8) :: IWaSrDm1Beg = 0.0d0
		real(8) :: IWaSrDm2Beg = 0.0d0

		! Per-compartment arrays
		real(8), allocatable :: DiPoCp(:)      !! Polygon diameter
		real(8), allocatable :: FrArMtrx(:)    !! Matrix area fraction
		real(8), allocatable :: VlMpDyCp(:)    !! Dynamic macropore volume
		real(8), allocatable :: VlMpStCp(:)    !! Static macropore volume
		real(8), allocatable :: VlMpStDm1(:)   !! Static volume domain 1
		real(8), allocatable :: VlMpStDm2(:)   !! Static volume domain 2
		real(8), allocatable :: QExcMpMtx(:)   !! Exchange flux
		real(8), allocatable :: SubsidCp(:)    !! Vertical subsidence
		real(8), allocatable :: PpDmCp(:,:)    !! Domain proportion
		real(8), allocatable :: dFdhMp(:)      !! Contribution to derivative
		real(8), allocatable :: iQOutDrRapCp(:)    !! Incremental rapid drainage per comp
		real(8), allocatable :: iQExcMtxDm1Cp(:)   !! Incremental exchange domain 1
		real(8), allocatable :: iQExcMtxDm2Cp(:)   !! Incremental exchange domain 2
		real(8), allocatable :: IAvFrMpWlWtDm1(:)  !! Avg wet wall fraction MB
		real(8), allocatable :: IAvFrMpWlWtDm2(:)  !! Avg wet wall fraction IC

		! Work arrays from SAVE (task persistence)
		integer, allocatable :: ICpBtDm(:)     !! Bottom compartment per domain
		integer, allocatable :: ICpTpWaSrDm(:) !! Top compartment water storage per domain
		real(8), allocatable :: ArMpTpDm(:)    !! Area at top per domain
		real(8), allocatable :: AwlCorFac(:)   !! Correction factor
		real(8), allocatable :: FrMpWalWet(:,:)    !! Wet macropore wall fraction
		real(8), allocatable :: KDCrRlRef(:)   !! Reference crack conductivity
		real(8), allocatable :: QExcMtxDmCp(:,:)   !! Exchange flux per domain/comp
		real(8), allocatable :: QInIntSatDmCp(:,:) !! Interflow in per domain/comp
		real(8), allocatable :: QInMtxSatDmCp(:,:) !! Matrix sat inflow per domain/comp
		real(8), allocatable :: QInTopLatDm(:)     !! Top lateral inflow per domain
		real(8), allocatable :: QInTopVrtDm(:)     !! Top vertical inflow per domain
		real(8), allocatable :: QOutDrRapCp(:)     !! Rapid drainage out per comp
		real(8), allocatable :: QOutMtxSatDmCp(:,:)    !! Matrix sat outflow per domain/comp
		real(8), allocatable :: QOutMtxUnsDmCp(:,:)    !! Matrix unsat outflow per domain/comp
		real(8), allocatable :: SorpDmCp(:,:)      !! Sorptivity per domain/comp
		real(8), allocatable :: ThtSrpRefDmCp(:,:) !! Reference theta for sorption
		real(8), allocatable :: TimAbsCumDmCp(:,:) !! Cumulative absorption time
		real(8), allocatable :: VlMpDm(:)          !! Macropore volume per domain
		real(8), allocatable :: VlMpDmCp(:,:)      !! Macropore volume per domain/comp
		real(8), allocatable :: WaSrMpDm(:)        !! Water storage per domain
		real(8), allocatable :: WaSrMpDmCp(:,:)    !! Water storage per domain/comp
		real(8), allocatable :: ZBtDm(:)           !! Bottom depth per domain
		real(8), allocatable :: ZWaLevDm(:)        !! Water level per domain
		logical, allocatable :: flDraTub(:)        !! Drain tube flag per level
		logical, allocatable :: FlEndSrpEvt(:,:)   !! End sorption event flag

		! State tracking
		real(8) :: WaSrMp = 0.0d0              !! Total water storage in macropores
		integer :: ICpBtPerZon = 0             !! Bottom compartment percolation zone
		integer :: ICpSatGWl = 0               !! Saturated compartment at GWL
		integer :: ICpSatPeGWl = 0             !! Saturated perched GWL compartment
		integer :: ICpTpPerZon = 0             !! Top compartment percolation zone
		integer :: ICpTpSatZon = 0             !! Top compartment saturated zone
		integer :: NnCrAr = 0                  !! Number of crack areas
		logical :: flBegin = .true.            !! Beginning of simulation flag

		! Domain configuration
		integer :: NumDm = 0                   !! Number of domains
		integer :: NumSbDm = 0                 !! Subdomains in IC
		integer :: IcTopMP = 0                 !! Top compartment with macropores
		integer :: NumLevRapDra = 0            !! Number of rapid drainage levels
		real(8) :: Z_Tp = 0.0d0                !! Top depth of macropores
		real(8) :: Z_St = 0.0d0                !! Bottom of static macropores
		real(8) :: Z_Ic = 0.0d0                !! Bottom of IC domain
		real(8) :: Z_Ah = 0.0d0                !! Bottom of A-horizon
		real(8) :: ArMpTp = 0.0d0              !! Area fraction at top of macropores
		real(8) :: ArMpSs = 0.0d0              !! Area fraction at soil surface
		real(8) :: KsatCovLay = 0.0d0          !! Saturated K of covering layer
		real(8) :: KsMpSs = 0.0d0              !! Vertical K of macropores at surface
		real(8) :: PpIcTpMp = 0.0d0            !! Proportion IC at top macropores
		real(8) :: dtold = 0.0d0               !! Previous timestep length

		! Groundwater tracking
		real(8) :: GWlFlCpZo = 0.0d0           !! GWL of full capillary zone
		integer :: NodGWlFlCpZo = 0
		real(8) :: ZDraBas = 0.0d0             !! Drainage basis level

		! Iteration control
		integer :: IDecMpRat = 0               !! Convergence iteration counter
		logical :: FlDecMpRat = .false.        !! Decrease macropore fluxes
		logical :: flInitDraBas = .false.      !! Initialize drainage basis

		! Flags
		logical :: flmacropore = .false.
	end type macropore_state_t

contains

	!> Initialize the macropore state.
	!!
	!! Allocates all dynamic arrays and resets fields to canonical defaults.
	!!
	!! @param[inout] state Macropore state container
	!! @param[in] numnod Number of soil nodes
	!! @param[in] maxdom Maximum number of macropore domains
	!! @param[in] maxdra Maximum number of drainage levels
	subroutine macropore_state_init(state, numnod, maxdom, maxdra)
		type(macropore_state_t), intent(inout) :: state
		integer, intent(in) :: numnod, maxdom, maxdra

		allocate(state%DiPoCp(numnod))
		allocate(state%FrArMtrx(numnod))
		allocate(state%VlMpDyCp(numnod))
		allocate(state%VlMpStCp(numnod))
		allocate(state%VlMpStDm1(numnod))
		allocate(state%VlMpStDm2(numnod))
		allocate(state%QExcMpMtx(numnod))
		allocate(state%SubsidCp(numnod))
		allocate(state%PpDmCp(maxdom, numnod))
		allocate(state%dFdhMp(numnod))
		allocate(state%iQOutDrRapCp(numnod))
		allocate(state%iQExcMtxDm1Cp(numnod))
		allocate(state%iQExcMtxDm2Cp(numnod))
		allocate(state%IAvFrMpWlWtDm1(numnod))
		allocate(state%IAvFrMpWlWtDm2(numnod))

		allocate(state%ICpBtDm(maxdom))
		allocate(state%ICpTpWaSrDm(maxdom))
		allocate(state%ArMpTpDm(maxdom))
		allocate(state%AwlCorFac(numnod))
		allocate(state%FrMpWalWet(maxdom, numnod))
		allocate(state%KDCrRlRef(maxdra))
		allocate(state%QExcMtxDmCp(maxdom, numnod))
		allocate(state%QInIntSatDmCp(maxdom, numnod))
		allocate(state%QInMtxSatDmCp(maxdom, numnod))
		allocate(state%QInTopLatDm(maxdom))
		allocate(state%QInTopVrtDm(maxdom))
		allocate(state%QOutDrRapCp(numnod))
		allocate(state%QOutMtxSatDmCp(maxdom, numnod))
		allocate(state%QOutMtxUnsDmCp(maxdom, numnod))
		allocate(state%SorpDmCp(maxdom, numnod))
		allocate(state%ThtSrpRefDmCp(maxdom, numnod))
		allocate(state%TimAbsCumDmCp(maxdom, numnod))
		allocate(state%VlMpDm(maxdom))
		allocate(state%VlMpDmCp(maxdom, numnod))
		allocate(state%WaSrMpDm(maxdom))
		allocate(state%WaSrMpDmCp(maxdom, numnod))
		allocate(state%ZBtDm(maxdom))
		allocate(state%ZWaLevDm(maxdom))
		allocate(state%flDraTub(maxdra))
		allocate(state%FlEndSrpEvt(maxdom, numnod))

		call macropore_state_reset_all(state)

		! Default matrix area fraction starts fully matrix-dominated
		state%FrArMtrx = 1.0d0

	end subroutine macropore_state_init

	!> Finalize the macropore state.
	!!
	!! Deallocates dynamic arrays and resets scalar fields.
	!!
	!! @param[inout] state Macropore state container
	subroutine macropore_state_finalize(state)
		type(macropore_state_t), intent(inout) :: state

		if (allocated(state%DiPoCp)) deallocate(state%DiPoCp)
		if (allocated(state%FrArMtrx)) deallocate(state%FrArMtrx)
		if (allocated(state%VlMpDyCp)) deallocate(state%VlMpDyCp)
		if (allocated(state%VlMpStCp)) deallocate(state%VlMpStCp)
		if (allocated(state%VlMpStDm1)) deallocate(state%VlMpStDm1)
		if (allocated(state%VlMpStDm2)) deallocate(state%VlMpStDm2)
		if (allocated(state%QExcMpMtx)) deallocate(state%QExcMpMtx)
		if (allocated(state%SubsidCp)) deallocate(state%SubsidCp)
		if (allocated(state%PpDmCp)) deallocate(state%PpDmCp)
		if (allocated(state%dFdhMp)) deallocate(state%dFdhMp)
		if (allocated(state%iQOutDrRapCp)) deallocate(state%iQOutDrRapCp)
		if (allocated(state%iQExcMtxDm1Cp)) deallocate(state%iQExcMtxDm1Cp)
		if (allocated(state%iQExcMtxDm2Cp)) deallocate(state%iQExcMtxDm2Cp)
		if (allocated(state%IAvFrMpWlWtDm1)) deallocate(state%IAvFrMpWlWtDm1)
		if (allocated(state%IAvFrMpWlWtDm2)) deallocate(state%IAvFrMpWlWtDm2)

		if (allocated(state%ICpBtDm)) deallocate(state%ICpBtDm)
		if (allocated(state%ICpTpWaSrDm)) deallocate(state%ICpTpWaSrDm)
		if (allocated(state%ArMpTpDm)) deallocate(state%ArMpTpDm)
		if (allocated(state%AwlCorFac)) deallocate(state%AwlCorFac)
		if (allocated(state%FrMpWalWet)) deallocate(state%FrMpWalWet)
		if (allocated(state%KDCrRlRef)) deallocate(state%KDCrRlRef)
		if (allocated(state%QExcMtxDmCp)) deallocate(state%QExcMtxDmCp)
		if (allocated(state%QInIntSatDmCp)) deallocate(state%QInIntSatDmCp)
		if (allocated(state%QInMtxSatDmCp)) deallocate(state%QInMtxSatDmCp)
		if (allocated(state%QInTopLatDm)) deallocate(state%QInTopLatDm)
		if (allocated(state%QInTopVrtDm)) deallocate(state%QInTopVrtDm)
		if (allocated(state%QOutDrRapCp)) deallocate(state%QOutDrRapCp)
		if (allocated(state%QOutMtxSatDmCp)) deallocate(state%QOutMtxSatDmCp)
		if (allocated(state%QOutMtxUnsDmCp)) deallocate(state%QOutMtxUnsDmCp)
		if (allocated(state%SorpDmCp)) deallocate(state%SorpDmCp)
		if (allocated(state%ThtSrpRefDmCp)) deallocate(state%ThtSrpRefDmCp)
		if (allocated(state%TimAbsCumDmCp)) deallocate(state%TimAbsCumDmCp)
		if (allocated(state%VlMpDm)) deallocate(state%VlMpDm)
		if (allocated(state%VlMpDmCp)) deallocate(state%VlMpDmCp)
		if (allocated(state%WaSrMpDm)) deallocate(state%WaSrMpDm)
		if (allocated(state%WaSrMpDmCp)) deallocate(state%WaSrMpDmCp)
		if (allocated(state%ZBtDm)) deallocate(state%ZBtDm)
		if (allocated(state%ZWaLevDm)) deallocate(state%ZWaLevDm)
		if (allocated(state%flDraTub)) deallocate(state%flDraTub)
		if (allocated(state%FlEndSrpEvt)) deallocate(state%FlEndSrpEvt)

		call macropore_state_reset_all(state)

	end subroutine macropore_state_finalize

	!> Reset cumulative macropore accumulators.
	!!
	!! @param[inout] state Macropore state container
	subroutine macropore_state_reset_cumulative(state)
		type(macropore_state_t), intent(inout) :: state

		state%cQMpLatSs = 0.0d0
		state%cQMpOutDrRap = 0.0d0
		state%cQMpInMtxSatDm1 = 0.0d0
		state%cQMpInMtxSatDm2 = 0.0d0
		state%cQMpOutMtxUnsDm1 = 0.0d0
		state%cQMpOutMtxUnsDm2 = 0.0d0
		state%cQMpInIntSatDm1 = 0.0d0
		state%cQMpInIntSatDm2 = 0.0d0
		state%cQMpInTopLatDm1 = 0.0d0
		state%cQMpInTopLatDm2 = 0.0d0
		state%cQMpInTopVrtDm1 = 0.0d0
		state%cQMpInTopVrtDm2 = 0.0d0
		state%cQMpOutMtxSatDm1 = 0.0d0
		state%cQMpOutMtxSatDm2 = 0.0d0

	end subroutine macropore_state_reset_cumulative

	!> Reset intermediate macropore accumulators.
	!!
	!! @param[inout] state Macropore state container
	subroutine macropore_state_reset_intermediate(state)
		type(macropore_state_t), intent(inout) :: state

		state%iQMpOutDrRap = 0.0d0
		state%iQInTopLatDm1 = 0.0d0
		state%iQInTopLatDm2 = 0.0d0
		state%iQInTopVrtDm1 = 0.0d0
		state%iQInTopVrtDm2 = 0.0d0
		state%IWaSrDm1Beg = 0.0d0
		state%IWaSrDm2Beg = 0.0d0

		if (allocated(state%iQOutDrRapCp)) state%iQOutDrRapCp = 0.0d0
		if (allocated(state%iQExcMtxDm1Cp)) state%iQExcMtxDm1Cp = 0.0d0
		if (allocated(state%iQExcMtxDm2Cp)) state%iQExcMtxDm2Cp = 0.0d0
		if (allocated(state%IAvFrMpWlWtDm1)) state%IAvFrMpWlWtDm1 = 0.0d0
		if (allocated(state%IAvFrMpWlWtDm2)) state%IAvFrMpWlWtDm2 = 0.0d0

	end subroutine macropore_state_reset_intermediate

	!> Reset all macropore state fields to canonical defaults.
	!!
	!! @param[inout] state Macropore state container
	subroutine macropore_state_reset_all(state)
		type(macropore_state_t), intent(inout) :: state

		state%VlMp = 0.0d0
		state%VlMpDm1 = 0.0d0
		state%VlMpDm2 = 0.0d0
		state%WaSrDm1 = 0.0d0
		state%WaSrDm2 = 0.0d0
		state%WaSrDm1Ini = 0.0d0
		state%WaSrDm2Ini = 0.0d0
		state%WaLevDm1 = 0.0d0
		state%QMaPo = 0.0d0
		state%QRapDra = 0.0d0
		state%QMpLatSs = 0.0d0
		state%QInTopLatDm1 = 0.0d0
		state%QInTopLatDm2 = 0.0d0
		state%QInTopVrtDm1 = 0.0d0
		state%QInTopVrtDm2 = 0.0d0

		if (allocated(state%DiPoCp)) state%DiPoCp = 0.0d0
		if (allocated(state%FrArMtrx)) state%FrArMtrx = 0.0d0
		if (allocated(state%VlMpDyCp)) state%VlMpDyCp = 0.0d0
		if (allocated(state%VlMpStCp)) state%VlMpStCp = 0.0d0
		if (allocated(state%VlMpStDm1)) state%VlMpStDm1 = 0.0d0
		if (allocated(state%VlMpStDm2)) state%VlMpStDm2 = 0.0d0
		if (allocated(state%QExcMpMtx)) state%QExcMpMtx = 0.0d0
		if (allocated(state%SubsidCp)) state%SubsidCp = 0.0d0
		if (allocated(state%PpDmCp)) state%PpDmCp = 0.0d0
		if (allocated(state%dFdhMp)) state%dFdhMp = 0.0d0
		if (allocated(state%ICpBtDm)) state%ICpBtDm = 0
		if (allocated(state%ICpTpWaSrDm)) state%ICpTpWaSrDm = 0
		if (allocated(state%ArMpTpDm)) state%ArMpTpDm = 0.0d0
		if (allocated(state%AwlCorFac)) state%AwlCorFac = 0.0d0
		if (allocated(state%FrMpWalWet)) state%FrMpWalWet = 0.0d0
		if (allocated(state%KDCrRlRef)) state%KDCrRlRef = 0.0d0
		if (allocated(state%QExcMtxDmCp)) state%QExcMtxDmCp = 0.0d0
		if (allocated(state%QInIntSatDmCp)) state%QInIntSatDmCp = 0.0d0
		if (allocated(state%QInMtxSatDmCp)) state%QInMtxSatDmCp = 0.0d0
		if (allocated(state%QInTopLatDm)) state%QInTopLatDm = 0.0d0
		if (allocated(state%QInTopVrtDm)) state%QInTopVrtDm = 0.0d0
		if (allocated(state%QOutDrRapCp)) state%QOutDrRapCp = 0.0d0
		if (allocated(state%QOutMtxSatDmCp)) state%QOutMtxSatDmCp = 0.0d0
		if (allocated(state%QOutMtxUnsDmCp)) state%QOutMtxUnsDmCp = 0.0d0
		if (allocated(state%SorpDmCp)) state%SorpDmCp = 0.0d0
		if (allocated(state%ThtSrpRefDmCp)) state%ThtSrpRefDmCp = 0.0d0
		if (allocated(state%TimAbsCumDmCp)) state%TimAbsCumDmCp = 0.0d0
		if (allocated(state%VlMpDm)) state%VlMpDm = 0.0d0
		if (allocated(state%VlMpDmCp)) state%VlMpDmCp = 0.0d0
		if (allocated(state%WaSrMpDm)) state%WaSrMpDm = 0.0d0
		if (allocated(state%WaSrMpDmCp)) state%WaSrMpDmCp = 0.0d0
		if (allocated(state%ZBtDm)) state%ZBtDm = 0.0d0
		if (allocated(state%ZWaLevDm)) state%ZWaLevDm = 0.0d0
		if (allocated(state%flDraTub)) state%flDraTub = .false.
		if (allocated(state%FlEndSrpEvt)) state%FlEndSrpEvt = .false.

		state%WaSrMp = 0.0d0
		state%ICpBtPerZon = 0
		state%ICpSatGWl = 0
		state%ICpSatPeGWl = 0
		state%ICpTpPerZon = 0
		state%ICpTpSatZon = 0
		state%NnCrAr = 0
		state%flBegin = .true.

		state%NumDm = 0
		state%NumSbDm = 0
		state%IcTopMP = 0
		state%NumLevRapDra = 0
		state%Z_Tp = 0.0d0
		state%Z_St = 0.0d0
		state%Z_Ic = 0.0d0
		state%Z_Ah = 0.0d0
		state%ArMpTp = 0.0d0
		state%ArMpSs = 0.0d0
		state%KsatCovLay = 0.0d0
		state%KsMpSs = 0.0d0
		state%PpIcTpMp = 0.0d0
		state%dtold = 0.0d0

		state%GWlFlCpZo = 0.0d0
		state%NodGWlFlCpZo = 0
		state%ZDraBas = 0.0d0

		state%IDecMpRat = 0
		state%FlDecMpRat = .false.
		state%flInitDraBas = .false.
		state%flmacropore = .false.

		call macropore_state_reset_cumulative(state)
		call macropore_state_reset_intermediate(state)

	end subroutine macropore_state_reset_all

end module macropore_state_mod
