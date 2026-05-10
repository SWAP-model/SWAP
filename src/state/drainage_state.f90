!> @file drainage_state.f90
!! Typed state record for the drainage subsystem. Owns the 6
!! drainage-flux variables that legacy SWAP held in the
!! `variables.f90` globals module.
!!
!! `qdra` and `qdrain` were temporarily classified into
!! surfacewater_state_t during the surface-water migration pilot
!! (ADR 0030) because of write-site overlap. The drainage migration
!! arc moves them here, where they architecturally belong.
!!
!! See ADR 0031 (state-migration drainage subsystem) and the
!! 2026-05-10 design spec.

module drainage_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: drainage_state_t

   type :: drainage_state_t
      ! Scalar fluxes
      real(real64) :: qdrd     = 0.0_real64    ! drain-direction sub-flux (cm/d)

      ! Per-level (Madr-sized) arrays — allocated by drainage_init
      real(real64), allocatable :: qdrain(:)       ! lateral drainage flux per level (cm/d)
      real(real64), allocatable :: drainl(:)       ! drain length per level (cm)
      real(real64), allocatable :: wetper(:)       ! wetted perimeter per level (cm)
      real(real64), allocatable :: ztopdislay(:)   ! top of discharge layer per level (cm)

      ! Per-level / per-compartment array (Madr × macp)
      real(real64), allocatable :: qdra(:,:)       ! lateral drainage flux per level/compartment (cm/d)
   end type drainage_state_t

end module drainage_state_mod
