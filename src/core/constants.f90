!> Module containing numerical tolerance and limit constants for SWAP
!!
!! This module defines small threshold values and limits used throughout
!! the SWAP model for numerical comparisons, convergence criteria, and
!! boundary conditions.
!!
!! @note These values are defined with double precision (real(8)) for
!!       consistency across all SWAP computations.
!!
!! @author Original SWAP development team
!! @version $Id: params.fi 277 2016-01-28 20:48:55Z kroes006 $
module swap_constants
  implicit none
  public

  !> Nearly negligible threshold value (1.0e-10)
  !! Used for checking near-zero values and numerical noise
  real(8), parameter :: NIHIL = 1.0d-10

  !> Tiny threshold value (1.0e-3)
  !! Used for small but non-negligible comparisons
  real(8), parameter :: TINY = 1.0d-3

  !> Small threshold value (1.0e-6)
  !! Used for convergence criteria and precision limits
  real(8), parameter :: SMALL = 1.0d-6

  !> Very large value (1.0e12)
  !! Used as upper bound or sentinel value
  real(8), parameter :: VLARGE = 1.0d12

end module swap_constants
