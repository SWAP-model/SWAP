! File VersionID:
!   $Id: tridag.f90 341 2017-09-29 18:12:25Z kroes006 $
! ----------------------------------------------------------------------

!> Module providing numerical solvers for linear systems
!! This module contains routines for solving tridiagonal and band diagonal
!! linear systems of equations.
module numericalsolvers_mod
   use error_mod, only: fatalerr_collected
  implicit none
  public :: tridag, bandec, banbks
contains
  
  !> Solves a tridiagonal linear system of equations
  !! 
  !! Solves for a vector U in the tridiagonal linear system A*U = R using
  !! the Thomas algorithm (tridiagonal matrix algorithm).
  !!
  !! @param[in] n Number of equations
  !! @param[in] a Lower diagonal coefficients (size n)
  !! @param[in] b Main diagonal coefficients (size n)
  !! @param[in] c Upper diagonal coefficients (size n)
  !! @param[in] r Right-hand side vector (size n)
  !! @param[out] u Solution vector (size n)
  !! @param[out] ierror Error code (0 = success, non-zero = error)
  !!
  !! @note
  !! Date: 15/9/99
  !! 
  !! References:
  !! Press, W.H., B.P. Flannery, S.A. Teukolsky & W.T. Vetterling, 1989.
  !! Numerical Recipes in FORTRAN. Cambridge University Press, New York.
  !! pp 40-41
  !! @endnote
  subroutine TRIDAG(N, A, B, C, R, U, ierror)
    use swap_array_dimensions, only: macp
    implicit none
    
    ! Input arguments
    integer, intent(in)    :: n
    real(8), intent(in)    :: a(n), b(n), c(n), r(n)
    
    ! Output arguments
    real(8), intent(out)   :: u(n)
    integer, intent(out)   :: ierror
    
    ! Parameters
    real(8), parameter     :: small = 0.3d-37
    
    ! Local variables
    integer                :: i
    real(8)                :: gamma(macp), beta
    character(len=200)     :: messag
    
    ierror = 0
    
    ! (1) If b(1)=0 then rewrite the equations as a set of order n-1
    ! to eliminate u(2)
    if (abs(b(1)) .lt. small) then
      messag = 'During the numerical solution the factor b(1)' // &
               ' became too small !'
      ierror = 1000
      ! call fatalerr ('tridag',messag)
    else
      ! (2) Decomposition and forward substitution
      beta = b(1)
      u(1) = r(1) / beta
      do i = 2, n
        gamma(i) = c(i-1) / beta
        beta = b(i) - a(i)*gamma(i)
        
        ! (2.1) If beta=0 then go to another algorithm including
        ! elimination with pivoting
        if (abs(beta) .lt. small) then
          messag = 'during the numerical solution the factor beta' // &
                   ' became too small !'
          ierror = 1000 + i
          ! call fatalerr ('tridag',messag)
          messag = messag  ! for Forcheck
          return
        else
          u(i) = (r(i) - a(i)*u(i-1)) / beta
        end if
      end do
      
      ! (3) Back substitution
      do i = n-1, 1, -1
        u(i) = u(i) - gamma(i+1)*u(i+1)
      end do
    end if
    
    return
  end subroutine tridag
  
  !> LU decomposition of a band diagonal matrix
  !!
  !! Constructs an LU decomposition of a rowwise permutation of a band diagonal
  !! matrix A with m1 subdiagonal rows and m2 superdiagonal rows. The matrix is
  !! compactly stored in array a(1:n,1:m1+m2+1). The upper triangular matrix
  !! replaces a, while the lower triangular matrix is returned in al(1:n,1:m1).
  !!
  !! @param[inout] a Band diagonal matrix (np,mp), replaced by upper triangular matrix
  !! @param[in] n Number of equations
  !! @param[in] m1 Number of subdiagonal rows
  !! @param[in] m2 Number of superdiagonal rows
  !! @param[in] np First dimension of array a
  !! @param[in] mp Second dimension of array a
  !! @param[out] al Lower triangular matrix (np,mpl)
  !! @param[in] mpl Second dimension of array al
  !! @param[out] indx Output vector recording row permutation (n)
  !! @param[out] d Output as +1 or -1 depending on number of row interchanges
  !!
  !! @note
  !! Reference: Numerical Recipes, Chapter 2.4
  !! 
  !! This routine is used in combination with banbks to solve band-diagonal
  !! sets of equations.
  !! @endnote
  subroutine bandec(a, n, m1, m2, np, mp, al, mpl, indx, d)
    implicit none
    
    ! Input/output arguments
    integer, intent(in)    :: m1, m2, mp, mpl, n, np
    integer, intent(out)   :: indx(n)
    real(8), intent(out)   :: d
    real(8), intent(inout) :: a(np,mp)
    real(8), intent(out)   :: al(np,mpl)
    
    ! Parameters
    real(8), parameter     :: TINY = 1.d-20
    
    ! Local variables
    integer                :: i, j, k, l, mm
    real(8)                :: dum
    character(len=200)     :: messag
    
    mm = m1 + m2 + 1
    ! Check array dimensions
    if (mm .gt. mp .or. m1 .gt. mpl .or. n .gt. np) then
      messag = 'bad args in bandec !'
      call fatalerr_collected('bandec', messag)
    end if
    
    ! Rearrange the storage a bit
    l = m1
    do i = 1, m1
      do j = m1 + 2 - i, mm
        a(i, j-l) = a(i, j)
      end do
      l = l - 1
      do j = mm - l, mm
        a(i, j) = 0.0d0
      end do
    end do
    
    d = 1.0d0
    l = m1
    
    ! For each row...
    do k = 1, n
      dum = a(k, 1)
      i = k
      if (l .lt. n) l = l + 1
      
      ! Find the pivot element
      do j = k + 1, l
        if (abs(a(j,1)) .gt. abs(dum)) then
          dum = a(j, 1)
          i = j
        end if
      end do
      
      indx(k) = i
      ! Matrix is algorithmically singular, but proceed anyway with TINY pivot
      ! (desirable in some applications)
      if (abs(dum) .lt. 1.0d-20) a(k, 1) = TINY
      
      ! Interchange rows
      if (i .ne. k) then
        d = -d
        do j = 1, mm
          dum = a(k, j)
          a(k, j) = a(i, j)
          a(i, j) = dum
        end do
      end if
      
      ! Do the elimination
      do i = k + 1, l
        dum = a(i, 1) / a(k, 1)
        al(k, i-k) = dum
        do j = 2, mm
          a(i, j-1) = a(i, j) - dum * a(k, j)
        end do
        a(i, mm) = 0.0d0
      end do
    end do
    
    return
  end subroutine bandec
  
  !> Backsubstitution for band diagonal linear systems
  !!
  !! Given the arrays a, al, and indx as returned from bandec, and given a
  !! right-hand side vector b(1:n), solves the band diagonal linear equations
  !! A · x = b. The solution vector x overwrites b(1:n). The other input arrays
  !! are not modified, and can be left in place for successive calls with
  !! different right-hand sides.
  !!
  !! @param[in] a Upper triangular matrix from bandec (np,mp)
  !! @param[in] n Number of equations
  !! @param[in] m1 Number of subdiagonal rows
  !! @param[in] m2 Number of superdiagonal rows
  !! @param[in] np First dimension of array a
  !! @param[in] mp Second dimension of array a
  !! @param[in] al Lower triangular matrix from bandec (np,mpl)
  !! @param[in] mpl Second dimension of array al
  !! @param[in] indx Row permutation vector from bandec (n)
  !! @param[inout] b Right-hand side vector on input, solution vector on output (n)
  !!
  !! @note
  !! Reference: Numerical Recipes, Chapter 2.4
  !! 
  !! This routine is used in combination with bandec to solve band-diagonal
  !! sets of equations.
  !! @endnote
  subroutine banbks(a, n, m1, m2, np, mp, al, mpl, indx, b)
    implicit none
    
    ! Input arguments
    integer, intent(in)    :: m1, m2, mp, mpl, n, np
    integer, intent(in)    :: indx(n)
    real(8), intent(in)    :: a(np,mp)
    real(8), intent(in)    :: al(np,mpl)
    
    ! Input/output arguments
    real(8), intent(inout) :: b(n)
    
    ! Local variables
    integer                :: i, k, l, mm
    real(8)                :: dum
    character(len=200)     :: messag
    
    mm = m1 + m2 + 1
    ! Check array dimensions
    if (mm .gt. mp .or. m1 .gt. mpl .or. n .gt. np) then
      messag = 'bad args in banbks !'
      call fatalerr_collected('banbks', messag)
    end if
    
    ! Forward substitution, unscrambling the permuted rows as we go
    l = m1
    do k = 1, n
      i = indx(k)
      if (i .ne. k) then
        dum = b(k)
        b(k) = b(i)
        b(i) = dum
      end if
      if (l .lt. n) l = l + 1
      do i = k + 1, l
        b(i) = b(i) - al(k, i-k) * b(k)
      end do
    end do
    
    ! Backsubstitution
    l = 1
    do i = n, 1, -1
      dum = b(i)
      do k = 2, l
        dum = dum - a(i, k) * b(k+i-1)
      end do
      b(i) = dum / a(i, 1)
      if (l .lt. mm) l = l + 1
    end do
    
    return
  end subroutine banbks

end module numericalsolvers_mod