! File VersionID:
!   $Id: sptabulated.f90 366 2018-01-10 11:12:43Z kroes006 $
! ----------------------------------------------------------------------
module doln
!   logical, parameter :: do_ln_trans = .false.
   logical, parameter :: do_ln_trans = .true.
end module doln

module doTSPACK
!   logical, parameter :: use_TSPACK = .false.
   logical, parameter :: use_TSPACK = .true.

! making use of TSPACK
   integer,          parameter               :: Mt = 1000
   logical,                            save  :: per, unifrm
   integer,                            save  :: NCD, IENDC, IER, NIT, IENDC_K
   integer,          dimension(Mt),    save  :: ICFLG
   integer,          parameter               :: LWK = 7*Mt
   double precision,                   save  :: BMAX, SM, SMTOL
   double precision, dimension(5,Mt),  save  :: B
   double precision, dimension(LWK),   save  :: WK2
   end module doTSPACK

module TSPACK

! contains the following routines/functions:
! public:   TSPBI    TSVAL1
!
! all:
!       1) Level 1 modules
!
!            These are divided into two groups.
!
!       a)  The following modules return knots (in the parametric
!           case), knot derivatives, tension factors, and, in the
!           case of smoothing, knot function values, which define
!           the fitting function (or functions in the parametric
!           case).  The naming convention should be evident from the
!           descriptions.
!        TSPSI    TSPSS    TSPSP    TSPBI    TSPBP
!
!       b)  The following modules return values, derivatives, or
!           integrals of the fitting function(s).
!        TSVAL1   TSVAL2   TSVAL3   TSINTL
!
!       2) Level 2 modules
!       a)  The following modules are called by TSPSP and TSPBP to
!           obtain a sequence of knots (parameter values) associated
!           with a parametric curve.  For some data sets, it might
!           be advantageous to replace these with routines that
!           implement an alternative method of parameterization.
!        ARCL2D   ARCL3D
!
!       b)  The following modules are called by the level 1, group (a)
!           modules to obtain knot derivatives (and values in the case
!           of SMCRV).
!        YPC1   YPC1P    YPC2  YPC2P    SMCRV
!
!       c)  The following modules are called by the level 1, group (a)
!           modules to obtain tension factors associated with knot
!           intervals.
!        SIGS   SIGBI    SIGBP
!
!       d)  The following functions are called by the level 1, group
!           (b) modules to obtain values and derivatives.  These pro-
!           vide a more convenient alternative to the level 1 routines
!           when a single value is needed.
!        HVAL   HPVAL    HPPVAL
!
!       3) Level 3 modules
!       a)  The following functions are called by SIGBI to compute
!           tension factors, and are convenient for obtaining an
!           optimal tension factor associated with a single interval.
!        SIG0   SIG1     SIG2
!
!       b)  The following modules are of general utility.
!        INTRVL SNHCSH   STORE
!
!       c)  The remaining modules are listed below.
!        B2TRI  B2TRIP   ENDSLP   YPCOEF
!

!               TSPACK:  Tension Spline Curve Fitting Package
!
!                              Robert J. Renka
!                                 05/27/91
!
!
!       I.  INTRODUCTION
!
!
!            The primary purpose of TSPACK is to construct a smooth
!       function which interpolates a discrete set of data points.
!       The function may be required to have either one or two con-
!       tinuous derivatives, and, in the C-2 case, several options
!       are provided for selecting end conditions.  If the accuracy
!       of the data does not warrant interpolation, a smoothing func-
!       tion (which does not pass through the data points) may be
!       constructed instead.  The fitting method is designed to avoid
!       extraneous inflection points (associated with rapidly varying
!       data values) and preserve local shape properties of the data
!       (monotonicity and convexity), or to satisfy the more general
!       constraints of bounds on function values or first derivatives.
!       The package also provides a parametric representation for con-
!       structing general planar curves and space curves.
!
!            The fitting function h(x) (or each component h(t) in the
!       case of a parametric curve) is defined locally, on each
!       interval associated with a pair of adjacent abscissae (knots),
!       by its values and first derivatives at the endpoints of the
!       interval, along with a nonnegative tension factor SIGMA
!       associated with the interval (h is a Hermite interpolatory
!       tension spline).  With SIGMA = 0, h is the cubic function
!       defined by the endpoint values and derivatives, and, as SIGMA
!       increases, h approaches the linear interpolant of the endpoint
!       values.  Since the linear interpolant preserves positivity,
!       monotonicity, and convexity of the data, h can be forced to
!       preserve these properties by choosing SIGMA sufficiently
!       large.  Also, since SIGMA varies with intervals, no more
!       tension than necessary is used in each interval, resulting in
!       a better fit and greater efficiency than is achieved with a
!       single constant tension factor.
!
!
!
!       II.  USAGE
!
!
!            TSPACK must be linked to a driver program which re-
!       serves storage, reads a data set, and calls the appropriate
!       subprograms selected from those described below in section
!       III.B.  Header comments in the software modules provide
!       details regarding the specification of input parameters and
!       the work space requirements.  It is recommended that curves
!       be plotted in order to assess their appropriateness for the
!       application.  This requires a user-supplied graphics package.
!
!
!
!       III.  SOFTWARE
!
!
!       A)  Code
!
!            The code is written in 1977 ANSI Standard Fortran.  All
!       variable and array names conform to the default typing con-
!       vention:  I-N for type INTEGER and A-H or O-Z for type REAL.
!       (There are no DOUBLE PRECISION variables.)  There are 32
!       modules, and they are ordered alphabetically in the package.
!       Each consists of the following sections:
!
!           1)  the module name and parameter list with spaces sepa-
!               rating the parameters into one to three subsets:
!               input parameters, I/O parameters, and output parame-
!               ters (in that order);
!           2)  type statements in which all parameters are typed
!               and arrays are dimensioned;
!           3)  a heading with the name of the package, identifica-
!               tion of the author, and date of the most recent
!               modification to the module;
!           4)  a description of the module's purpose and other rel-
!               evant information for the user;
!           5)  input parameter descriptions and output parameter
!               descriptions in the same order as the parameter
!               list;
!           6)  a list of other modules required (called either
!               directly or indirectly);
!           7)  a list of intrinsic functions called, if any; and
!           8)  the code, including comments.
!
!            Note that it is assumed that floating point underflow
!       results in assignment of the value zero.  If not the default,
!       this may be specified as either a compiler option or an
!       operating system option.  Also, overflow is avoided by re-
!       stricting arguments to the exponential function EXP to have
!       value at most SBIG=85.  SBIG, which appears in DATA statements
!       in the evaluation functions, HVAL, HPVAL, HPPVAL, and TSINTL,
!       must be decreased if it is necessary to accomodate a floating
!       point number system with fewer than 8 bits in the exponent.
!       No other system dependencies are present in the code.
!
!            The modules that solve nonlinear equations, SIGS, SIGBP,
!       SIG0, SIG1, SIG2, and SMCRV, include diagnostic print capabi-
!       lity which allows the iteration to be traced.  This can be
!       enabled by altering logical unit number LUN in a DATA state-
!       ment in the relevant module.
!
!
!       B)  Module Descriptions
!
!            The software modules are divided into three categories,
!       referred to as level 1, level 2, and level 3, corresponding
!       to the hierarchy of calling sequences:  level 1 modules call
!       level 2 modules which call level 3 modules.  For most ap-
!       plications, the driver need only call two level 1 modules --
!       one from each of groups (a) and (b).  However, additional
!       control over various options can be obtained by directly
!       calling level 2 modules.  Also, additional fitting methods,
!       such as parametric smoothing, can be obtained by calling
!       level 2 modules.  Note that, in the case of smoothing or C-2
!       interpolation with automatically selected tension, the use
!       of level 2 modules requires that an iteration be placed around
!       the computation of knot derivatives and tension factors.
!
!
!
!       1) Level 1 modules
!
!            These are divided into two groups.
!
!       a)  The following modules return knots (in the parametric
!           case), knot derivatives, tension factors, and, in the
!           case of smoothing, knot function values, which define
!           the fitting function (or functions in the parametric
!           case).  The naming convention should be evident from the
!           descriptions.
!
!
!       TSPSI   Subroutine which constructs a shape-preserving or
!                 unconstrained interpolatory function.  Refer to
!                 TSVAL1.
!
!       TSPSS   Subroutine which constructs a shape-preserving or
!                 unconstrained smoothing spline.  Refer to TSVAL1.
!
!       TSPSP   Subroutine which constructs a shape-preserving or
!                 unconstrained planar curve or space curve.  Refer
!                 to TSVAL2 and TSVAL3.
!
!       TSPBI   Subroutine which constructs a bounds-constrained
!                 interpolatory function.  The constraints are defined
!                 by a user-supplied array containing upper and lower
!                 bounds on function values and first derivatives,
!                 along with required signs for the second derivative,
!                 for each interval.  Refer to TSVAL1.
!
!       TSPBP   Subroutine which constructs a bounds-constrained
!                 parametric planar curve.  The constraints are de-
!                 fined by user-supplied arrays containing upper and
!                 lower bounds on the signed perpendicular distance
!                 between the smooth curve segment and the polygonal
!                 line segment associated with each knot interval.
!                 The bounds might, for example, be chosen to avoid
!                 intersections between smooth contour curves.  Refer
!                 to TSVAL2.
!
!
!       b)  The following modules return values, derivatives, or
!           integrals of the fitting function(s).
!
!
!       TSVAL1  Subroutine which evaluates a Hermite interpolatory
!                 tension spline or its first or second derivative at
!                 a user-specified set of points.  Note that a smooth-
!                 ing curve constructed by TSPSS is the interpolant of
!                 the computed knot function values.  The evaluation
!                 points need not lie in the interval defined by the
!                 knots, but care must be exercised in assessing the
!                 accuracy of extrapolation.
!
!       TSVAL2  Subroutine which returns values or derivatives of a
!                 pair of Hermite interpolatory tension splines which
!                 form the components of a parametric planar curve.
!                 The output values may be used to construct unit tan-
!                 gent vectors, curvature vectors, etc.
!
!       TSVAL3  Subroutine which returns values or derivatives of
!                 three Hermite interpolatory tension splines which
!                 form the components of a parametric space curve.
!
!       TSINTL  Function which returns the integral over a specified
!                 domain of a Hermite interpolatory tension spline.
!                 This provides an effective means of quadrature for
!                 a function defined only by a discrete set of values.
!
!
!
!       2) Level 2 modules
!
!            These are divided into four groups.
!
!       a)  The following modules are called by TSPSP and TSPBP to
!           obtain a sequence of knots (parameter values) associated
!           with a parametric curve.  For some data sets, it might
!           be advantageous to replace these with routines that
!           implement an alternative method of parameterization.
!
!
!       ARCL2D  Subroutine which computes the sequence of cumulative
!                 arc lengths associated with a sequence of points in
!                 the plane.
!
!       ARCL3D  Subroutine which computes the sequence of cumulative
!                 arc lengths associated with a sequence of points in
!                 3-space.
!
!
!       b)  The following modules are called by the level 1, group (a)
!           modules to obtain knot derivatives (and values in the case
!           of SMCRV).
!
!
!    YPC1       Subroutine which employs a monotonicity-constrained
!                quadratic interpolation method to compute locally
!                defined derivative estimates, resulting in a C-1
!                fit.
!
!       YPC1P    Subroutine similar to YPC1 for the case of periodic
!                 end conditions.  In the case of a parametric curve
!                 fit, periodic end conditions are necessary to ob-
!                 tain a closed curve.
!
!       YPC2     Subroutine which determines a set of knot-derivative
!                 estimates which result in a tension spline with two
!                 continuous derivatives and satisfying specified
!                 end conditions.
!
!       YPC2P    Subroutine similar to YPC2 for the case of periodic
!                 end conditions.
!
!       SMCRV    Subroutine which, given a sequence of abscissae with
!                 associated data values and tension factors, returns
!                 a set of function values and first derivative values
!                 defining a twice-continuously differentiable tension
!                spline which smoothes the data and satisfies either
!                natural or periodic end conditions.
!
!
!       c)  The following modules are called by the level 1, group (a)
!           modules to obtain tension factors associated with knot
!           intervals.
!
!
!       SIGS    Subroutine which, given a sequence of abscissae,
!                 function values, and first derivative values,
!                 determines the set of minimum tension factors for
!                 which the Hermite interpolatory tension spline
!                 preserves local shape properties (monotonicity
!                 and convexity) of the data.  SIGS is called by
!                 TSPSI, TSPSS, and TSPSP.
!
!       SIGBI   Subroutine which, given a sequence of abscissae,
!                 function values, and first derivative values,
!                 determines the set of minimum tension factors for
!                 which the Hermite interpolatory tension spline
!                 satisfies specified bounds constraints.  SIGBI is
!                 called by TSPBI.
!
!       SIGBP   Subroutine which, given a sequence of points in the
!                 plane with associated derivative vectors, determines
!                 the set of minimum tension factors for which the
!                 parametric planar tension spline curve defined by
!                 the data satisfies specified bounds on the signed
!                 orthogonal distance between the parametric curve and
!                 the polygonal curve defined by the points.  SIGBP is
!                 called by TSPBP.
!
!
!       d)  The following functions are called by the level 1, group
!           (b) modules to obtain values and derivatives.  These pro-
!           vide a more convenient alternative to the level 1 routines
!           when a single value is needed.
!
!
!       HVAL     Function which evaluates a Hermite interpolatory ten-
!                 sion spline at a specified point.
!
!       HPVAL    Function which evaluates the first derivative of a
!                 Hermite interpolatory tension spline at a specified
!                 point.
!
!       HPPVAL  Function which evaluates the second derivative of a
!                 Hermite interpolatory tension spline at a specified
!                 point.
!
!
!
!       3) Level 3 modules
!
!            These are divided into three groups.
!
!       a)  The following functions are called by SIGBI to compute
!           tension factors, and are convenient for obtaining an
!           optimal tension factor associated with a single interval.
!
!
!       SIG0     Function which, given a pair of abscissae, along with
!                 associated endpoint values and derivatives, deter-
!                 mines the smallest tension factor for which the
!                 corresponding Hermite interpolatory tension spline
!                 satisfies a specified bound on function values in
!                 the interval.
!
!       SIG1     Function which, given a pair of abscissae, along with
!                 associated endpoint data, determines the smallest
!                 tension factor for which the corresponding Hermite
!                 interpolatory tension spline satisfies a specified
!                 bound on first derivative values in the interval.
!
!       SIG2     Function which, given a pair of abscissae, along with
!                 associated endpoint data, determines the smallest
!                 tension factor for which the corresponding Hermite
!                 interpolatory tension spline preserves convexity
!                 (or concavity) of the data.
!
!
!       b)  The following modules are of general utility.
!
!
!       INTRVL Function which, given an increasing sequence of ab-
!                 scissae, returns the index of an interval containing
!                a specified point.  INTRVL is called by the evalua-
!                tion functions TSINTL, HVAL, HPVAL, and HPPVAL.
!
!       SNHCSH  Subroutine called by several modules to compute
!                 accurate approximations to the modified hyperbolic
!                 functions which form a basis for exponential ten-
!                 sion splines.
!
!       STORE   Function used by SIGBP, SIGS, SIG0, SIG1, and SIG2 in
!                 computing the machine precision.  STORE forces a
!                 value to be stored in main memory so that the pre-
!                 cision of floating point numbers in memory locations
!                 rather than registers is computed.
!
!
!       c)  The remaining modules are listed below.
!
!
!       B2TRI   Subroutine called by SMCRV to solve the symmetric
!                 positive-definite block tridiagonal linear system
!                 associated with the gradient of the quadratic
!                 functional whose minimum corresponds to a smooth-
!                 ing curve with nonperiodic end conditions.
!
!       B2TRIP  Subroutine similar to B2TRI for periodic end
!                 conditions.
!
!       ENDSLP  Function which returns the derivative at X1 of a
!                 tension spline h(x) which interpolates three
!                 specified data points and has third derivative
!                 equal to zero at X1.  ENDSLP is called by YPC2
!                 when this choice of end conditions is selected
!                 by an input parameter.
!
!       YPCOEF  Subroutine called by SMCRV, YPC2, and YPC2P to com-
!                 pute coefficients defining the linear system.
!
!
!
!       IV.  REFERENCE
!
!
!       For the theoretical background, consult the following:
!
!         RENKA, R. J.  Interpolatory tension splines with automatic
!         selection of tension factors. SIAM J. Sci. Stat. Comput. 8
!         (1987), pp. 393-415.
private
public :: TSPBI, TSVAL1, my_HVAL, my_HPVAL
contains

      DOUBLE PRECISION FUNCTION ENDSLP (X1,X2,X3,Y1,Y2,Y3,              &
     &                                  SIGMA)
      DOUBLE PRECISION X1, X2, X3, Y1, Y2, Y3, SIGMA
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/17/96
!C
!C   Given data values associated with a strictly increasing
!C or decreasing sequence of three abscissae X1, X2, and X3,
!C this function returns a derivative estimate at X1 based on
!C the tension spline H(x) which interpolates the data points
!C and has third derivative equal to zero at X1.  Letting S1
!C denote the slope defined by the first two points, the est-
!C mate is obtained by constraining the derivative of H at X1
!C so that it has the same sign as S1 and its magnitude is
!C at most 3*abs(S1).  If SIGMA = 0, H(x) is quadratic and
!C the derivative estimate is identical to the value computed
!C by Subroutine YPC1 at the first point (or the last point
!C if the abscissae are decreasing).
!C
!C On input:
!C
!C       X1,X2,X3 = Abscissae satisfying either X1 < X2 < X3
!C                  or X1 > X2 > X3.
!C
!C       Y1,Y2,Y3 = Data values associated with the abscis-
!C                  sae.  H(X1) = Y1, H(X2) = Y2, and H(X3)
!C                  = Y3.
!C
!C       SIGMA = Tension factor associated with H in inter-
!C               val (X1,X2) or (X2,X1).
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       ENDSLP = (Constrained) derivative of H at X1, or
!C                zero if the abscissae are not strictly
!C                monotonic.
!C
!C Module required by ENDSLP:  SNHCSH
!C
!C Intrinsic functions called by ENDSLP:  ABS, EXP, MAX, MIN
!C
!C***********************************************************
!C
      DOUBLE PRECISION COSHM1, COSHMS, DUMMY, DX1, DXS, S1,             &
     &                 SIG1, SIGS, T
!C
      DX1 = X2 - X1
      DXS = X3 - X1
      IF (DX1*(DXS-DX1) .LE. 0.D0) GO TO 2
      SIG1 = DABS(SIGMA)
      IF (SIG1 .LT. 1.D-9) THEN
!C
!C SIGMA = 0:  H is the quadratic interpolant.
!C
        T = (DX1/DXS)**2
        GO TO 1
      ENDIF
      SIGS = SIG1*DXS/DX1
      IF (SIGS .LE. .5D0) THEN
!C
!C 0 < SIG1 < SIGS .LE. .5:  compute approximations to
!C   COSHM1 = COSH(SIG1)-1 and COSHMS = COSH(SIGS)-1.
!C
        CALL SNHCSH (SIG1, DUMMY,COSHM1,DUMMY)
        CALL SNHCSH (SIGS, DUMMY,COSHMS,DUMMY)
        T = COSHM1/COSHMS
      ELSE
!C
!C SIGS > .5:  compute T = COSHM1/COSHMS.
!C
        T = EXP(SIG1-SIGS)*((1.D0-EXP(-SIG1))/                          &
     &                      (1.D0-EXP(-SIGS)))**2
      ENDIF
!C
!C The derivative of H at X1 is
!C   T = ((Y3-Y1)*COSHM1-(Y2-Y1)*COSHMS)/
!C       (DXS*COSHM1-DX1*COSHMS).
!C
!C ENDSLP = T unless T*S1 < 0 or abs(T) > 3*abs(S1).
!C
    1 T = ((Y3-Y1)*T-Y2+Y1)/(DXS*T-DX1)
      S1 = (Y2-Y1)/DX1
      IF (S1 .GE. 0.D0) THEN
        ENDSLP = MIN(MAX(0.D0,T), 3.D0*S1)
      ELSE
        ENDSLP = MAX(MIN(0.D0,T), 3.D0*S1)
      ENDIF
      RETURN
!C
!C Error in the abscissae.
!C
    2 ENDSLP = 0.D0
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION HPPVAL (T,N,X,Y,YP,                     &
     &                                  SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/17/96
!C
!C   This function evaluates the second derivative HPP of a
!C Hermite interpolatory tension spline H at a point T.
!C
!C On input:
!C
!C       T = Point at which HPP is to be evaluated.  Extrap-
!C           olation is performed if T < X(1) or T > X(N).
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing the abscissae.
!C           These must be in strictly increasing order:
!C           X(I) < X(I+1) for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values.
!C           H(X(I)) = Y(I) for I = 1,...,N.
!C
!C       YP = Array of length N containing first deriva-
!C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
!C            HP denotes the derivative of H.
!C
!C       SIGMA = Array of length N-1 containing tension fac-
!C               tors whose absolute values determine the
!C               balance between cubic and linear in each
!C               interval.  SIGMA(I) is associated with int-
!C               erval (I,I+1) for I = 1,...,N-1.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered and
!C                     X(1) .LE. T .LE. X(N).
!C             IER = 1 if no errors were encountered and
!C                     extrapolation was necessary.
!C             IER = -1 if N < 2.
!C             IER = -2 if the abscissae are not in strictly
!C                      increasing order.  (This error will
!C                      not necessarily be detected.)
!C
!C       HPPVAL = Second derivative value HPP(T), or zero if
!C                IER < 0.
!C
!C Modules required by HPPVAL:  INTRVL, SNHCSH
!C
!C Intrinsic functions called by HPPVAL:  ABS, EXP
!C
!C***********************************************************
!C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, COSH2, D1, D2,             &
     &                 DUMMY, DX, E, E1, E2, EMS, S, SB1,               &
     &                 SB2, SBIG, SIG, SINH2, SM, SM2, TM
! in a module these need not to be declared here      INTEGER INTRVL
!C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
!C
!C Find the index of the left end of an interval containing
!C   T.  If T < X(1) or T > X(N), extrapolation is performed
!C   using the leftmost or rightmost interval.
!C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
!C
!C Compute interval width DX, local coordinates B1 and B2,
!C   and second differences D1 and D2.
!C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      B1 = (X(IP1) - T)/DX
      B2 = 1.D0 - B1
      S = (Y(IP1)-Y(I))/DX
      D1 = S - YP(I)
      D2 = YP(IP1) - S
      SIG = DABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
!C
!C SIG = 0:  H is the Hermite cubic interpolant.
!C
        HPPVAL = (D1 + D2 + 3.D0*(B2-B1)*(D2-D1))/DX
      ELSEIF (SIG .LE. .5D0) THEN
!C
!C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!C   cancellation error in the hyperbolic functions.
!C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH2 = SM2 + SB2
        COSH2 = CM2 + 1.D0
        E = SIG*SM - CMM - CMM
        HPPVAL = SIG*((CM*SINH2-SM*COSH2)*(D1+D2) +                     &
     &              SIG*(CM*COSH2-(SM+SIG)*SINH2)*D1)/(DX*E)
      ELSE
!C
!C SIG > .5:  use negative exponentials in order to avoid
!C   overflow.  Note that EMS = EXP(-SIG).  In the case of
!C   extrapolation (negative B1 or B2), H is approximated by
!C   a linear function if -SIG*B1 or -SIG*B2 is large.
!C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HPPVAL = 0.D0
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          E = TM*(SIG*(1.D0+EMS) - TM - TM)
          HPPVAL = SIG*(SIG*((E1*EMS+E2)*D1+(E1+E2*EMS)*D2)-            &
     &                       TM*(E1+E2)*(D1+D2))/(DX*E)
        ENDIF
      ENDIF
      RETURN
!C
!C N is outside its valid range.
!C
    1 HPPVAL = 0.D0
      IER = -1
      RETURN
!C
!C X(I) .GE. X(I+1).
!C
    2 HPPVAL = 0.D0
      IER = -2
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION HPVAL (T,N,X,Y,YP,                      &
     &                                 SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/17/96
!C
!C   This function evaluates the first derivative HP of a
!C Hermite interpolatory tension spline H at a point T.
!C
!C On input:
!C
!C       T = Point at which HP is to be evaluated.  Extrapo-
!C           lation is performed if T < X(1) or T > X(N).
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing the abscissae.
!C           These must be in strictly increasing order:
!C           X(I) < X(I+1) for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values.
!C           H(X(I)) = Y(I) for I = 1,...,N.
!C
!C       YP = Array of length N containing first deriva-
!C            tives.  HP(X(I)) = YP(I) for I = 1,...,N.
!C
!C       SIGMA = Array of length N-1 containing tension fac-
!C               tors whose absolute values determine the
!C               balance between cubic and linear in each
!C               interval.  SIGMA(I) is associated with int-
!C               erval (I,I+1) for I = 1,...,N-1.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0  if no errors were encountered and
!C                      X(1) .LE. T .LE. X(N).
!C             IER = 1  if no errors were encountered and
!C                      extrapolation was necessary.
!C             IER = -1 if N < 2.
!C             IER = -2 if the abscissae are not in strictly
!C                      increasing order.  (This error will
!C                      not necessarily be detected.)
!C
!C       HPVAL = Derivative value HP(T), or zero if IER < 0.
!C
!C Modules required by HPVAL:  INTRVL, SNHCSH
!C
!C Intrinsic functions called by HPVAL:  ABS, EXP
!C
!C***********************************************************
!C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,             &
     &                 DX, E, E1, E2, EMS, S, S1, SB1, SB2,             &
     &                 SBIG, SIG, SINH2, SM, SM2, TM
! in a module these need not to be declared here      INTEGER INTRVL
!C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
!C
!C Find the index of the left end of an interval containing
!C   T.  If T < X(1) or T > X(N), extrapolation is performed
!C   using the leftmost or rightmost interval.
!C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
!C
!C Compute interval width DX, local coordinates B1 and B2,
!C   and second differences D1 and D2.
!C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      B1 = (X(IP1) - T)/DX
      B2 = 1.D0 - B1
      S1 = YP(I)
      S = (Y(IP1)-Y(I))/DX
      D1 = S - S1
      D2 = YP(IP1) - S
      SIG = DABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
!C
!C SIG = 0:  H is the Hermite cubic interpolant.
!C
        HPVAL = S1 + B2*(D1 + D2 - 3.D0*B1*(D2-D1))
      ELSEIF (SIG .LE. .5D0) THEN
!C
!C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!C   cancellation error in the hyperbolic functions.
!C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH2 = SM2 + SB2
        E = SIG*SM - CMM - CMM
        HPVAL = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) +                       &
     &                SIG*(CM*SINH2-(SM+SIG)*CM2)*D1)/E
      ELSE
!C
!C SIG > .5:  use negative exponentials in order to avoid
!C   overflow.  Note that EMS = EXP(-SIG).  In the case of
!C   extrapolation (negative B1 or B2), H is approximated by
!C   a linear function if -SIG*B1 or -SIG*B2 is large.
!C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HPVAL = S
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          E = TM*(SIG*(1.D0+EMS) - TM - TM)
          HPVAL = S + (TM*((E2-E1)*(D1+D2) + TM*(D1-D2)) +              &
     &            SIG*((E1*EMS-E2)*D1 + (E1-E2*EMS)*D2))/E
        ENDIF
      ENDIF
      RETURN
!C
!C N is outside its valid range.
!C
    1 HPVAL = 0.D0
      IER = -1
      RETURN
!C
!C X(I) .GE. X(I+1).
!C
    2 HPVAL = 0.D0
      IER = -2
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION HVAL (T,N,X,Y,YP,SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu

!C                                                   11/17/96
!C
!C   This function evaluates a Hermite interpolatory tension
!C spline H at a point T.  Note that a large value of SIGMA
!C may cause underflow.  The result is assumed to be zero.
!C
!C   Given arrays X, Y, YP, and SIGMA of length NN, if T is
!C known to lie in the interval (X(I),X(J)) for some I < J,
!C a gain in efficiency can be achieved by calling this
!C function with N = J+1-I (rather than NN) and the I-th
!C components of the arrays (rather than the first) as par-
!C ameters.
!C
!C On input:
!C
!C       T = Point at which H is to be evaluated.  Extrapo-
!C           lation is performed if T < X(1) or T > X(N).
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing the abscissae.
!C           These must be in strictly increasing order:
!C           X(I) < X(I+1) for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values.
!C           H(X(I)) = Y(I) for I = 1,...,N.
!C
!C       YP = Array of length N containing first deriva-
!C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
!C            HP denotes the derivative of H.
!C
!C       SIGMA = Array of length N-1 containing tension fac-
!C               tors whose absolute values determine the
!C               balance between cubic and linear in each
!C               interval.  SIGMA(I) is associated with int-
!C               erval (I,I+1) for I = 1,...,N-1.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0  if no errors were encountered and
!C                      X(1) .LE. T .LE. X(N).
!C             IER = 1  if no errors were encountered and
!C                      extrapolation was necessary.
!C             IER = -1 if N < 2.
!C             IER = -2 if the abscissae are not in strictly
!C                      increasing order.  (This error will
!C                      not necessarily be detected.)
!C
!C       HVAL = Function value H(T), or zero if IER < 0.
!C
!C Modules required by HVAL:  INTRVL, SNHCSH
!C
!C Intrinsic functions called by HVAL:  ABS, EXP
!C
!C***********************************************************
!C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,             &
     &                 DX, E, E1, E2, EMS, S, S1, SB1, SB2,             &
     &                 SBIG, SIG, SM, SM2, TM, TP, TS, U, Y1
! in a module these need not to be declared here      INTEGER INTRVL
!C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
!C
!C Find the index of the left end of an interval containing
!C   T.  If T < X(1) or T > X(N), extrapolation is performed
!C   using the leftmost or rightmost interval.
!C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
!C
!C Compute interval width DX, local coordinates B1 and B2,
!C   and second differences D1 and D2.
!C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      U = T - X(I)
      B2 = U/DX
      B1 = 1.D0 - B2
      Y1 = Y(I)
      S1 = YP(I)
      S = (Y(IP1)-Y1)/DX
      D1 = S - S1
      D2 = YP(IP1) - S
      SIG = DABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
!C
!C SIG = 0:  H is the Hermite cubic interpolant.
!C
        HVAL = Y1 + U*(S1 + B2*(D1 + B1*(D1-D2)))
      ELSEIF (SIG .LE. .5D0) THEN
!C
!C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!C   cancellation error in the hyperbolic functions.
!C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HVAL = Y1 + S1*U + DX*((CM*SM2-SM*CM2)*(D1+D2) +                &
     &                         SIG*(CM*CM2-(SM+SIG)*SM2)*D1)            &
     &                         /(SIG*E)
      ELSE
!C
!C SIG > .5:  use negative exponentials in order to avoid
!C   overflow.  Note that EMS = EXP(-SIG).  In the case of
!C   extrapolation (negative B1 or B2), H is approximated by
!C   a linear function if -SIG*B1 or -SIG*B2 is large.
!C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HVAL = Y1 + S*U
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TS = TM*TM
          TP = 1.D0 + EMS
          E = TM*(SIG*TP - TM - TM)
          HVAL = Y1 + S*U + DX*(TM*(TP-E1-E2)*(D1+D2) + SIG*            &
     &                         ((E2+EMS*(E1-2.D0)-B1*TS)*D1+            &
     &                        (E1+EMS*(E2-2.D0)-B2*TS)*D2))/            &
     &                        (SIG*E)
        ENDIF
      ENDIF
      RETURN
!C
!C N is outside its valid range.
!C
    1 HVAL = 0.D0
      IER = -1
      RETURN
!C
!C X(I) .GE. X(I+1).
!C
    2 HVAL = 0.D0
      IER = -2
      RETURN
      END FUNCTION

      INTEGER FUNCTION INTRVL (T,N,X)
      INTEGER N
      DOUBLE PRECISION T, X(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This function returns the index of the left end of an
!C interval (defined by an increasing sequence X) which
!C contains the value T.  The method consists of first test-
!C ing the interval returned by a previous call, if any, and
!C then using a binary search if necessary.
!C
!C On input:
!C
!C       T = Point to be located.
!C
!C       N = Length of X.  N .GE. 2.
!C
!C       X = Array of length N assumed (without a test) to
!C           contain a strictly increasing sequence of
!C           values.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       INTRVL = Index I defined as follows:
!C
!C                  I = 1    if  T .LT. X(2) or N .LE. 2,
!C                  I = N-1  if  T .GE. X(N-1), and
!C                  X(I) .LE. T .LT. X(I+1) otherwise.
!C
!C Modules required by INTRVL:  None
!C
!C***********************************************************
!C
      INTEGER IH, IL, K
      DOUBLE PRECISION TT
!C
      SAVE IL
      DATA IL/1/
      TT = T
      IF (IL .GE. 1  .AND.  IL .LT. N) THEN
        IF (X(IL) .LE. TT  .AND.  TT .LT. X(IL+1)) GO TO 2
      ENDIF
!C
!C Initialize low and high indexes.
!C
      IL = 1
      IH = N
!C
!C Binary search:
!C
    1 IF (IH .LE. IL+1) GO TO 2
        K = (IL+IH)/2
        IF (TT .LT. X(K)) THEN
          IH = K
        ELSE
          IL = K
        ENDIF
        GO TO 1
!C
!C X(IL) .LE. T .LT. X(IL+1)  or  (T .LT. X(1) and IL=1)
!C                            or  (T .GE. X(N) and IL=N-1)
!C
    2 INTRVL = IL
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION SIG0 (X1,X2,Y1,Y2,Y1P,Y2P,              &
     &                                IFL,HBND,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, HBND, TOL
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/18/96
!C
!C   Given a pair of abscissae with associated ordinates and
!C slopes, this function determines the smallest (nonnega-
!C tive) tension factor SIGMA such that the Hermite interpo-
!C latory tension spline H(x), defined by SIGMA and the data,
!C is bounded (either above or below) by HBND for all x in
!C (X1,X2).
!C
!C On input:
!C
!C       X1,X2 = Abscissae.  X1 < X2.
!C
!C       Y1,Y2 = Values of H at X1 and X2.
!C
!C       Y1P,Y2P = Derivative values of H at X1 and X2.
!C
!C       IFL = Option indicator:
!C             IFL = -1 if HBND is a lower bound on H.
!C             IFL = 1 if HBND is an upper bound on H.
!C
!C       HBND = Bound on H.  If IFL = -1, HBND .LE. min(Y1,
!C              Y2).  If IFL = 1, HBND .GE. max(Y1,Y2).
!C
!C       TOL = Tolerance whose magnitude determines how close
!C             SIGMA is to its optimal value when nonzero
!C             finite tension is necessary and sufficient to
!C             satisfy the constraint.  For a lower bound,
!C             SIGMA is chosen so that HBND .LE. HMIN .LE.
!C             HBND + abs(TOL), where HMIN is the minimum
!C             value of H in the interval, and for an upper
!C             bound, the maximum of H satisfies HBND -
!C             abs(TOL) .LE. HMAX .LE. HBND.  Thus, the con-
!C             straint is satisfied but possibly with more
!C             tension than necessary.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered and the
!C                     constraint can be satisfied with fin-
!C                     ite tension.
!C             IER = 1 if no errors were encountered but in-
!C                     finite tension is required to satisfy
!C                     the constraint (e.g., IFL = -1, HBND
!C                     = Y1, and Y1P < 0.).
!C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
!C             IER = -2 if HBND is outside its valid range
!C                      on input.
!C
!C       SIG0 = Minimum tension factor defined above unless
!C              IER < 0, in which case SIG0 = -1.  If IER =
!C              1, SIG0 = 85, resulting in an approximation
!C              to the linear interpolant of the endpoint
!C              values.  Note, however, that SIG0 may be
!C              larger than 85 if IER = 0.
!C
!C Modules required by SIG0:  SNHCSH, STORE
!C
!C Intrinsic functions called by SIG0:  ABS, DBLE, EXP, LOG,
!C                                        MAX, MIN, SIGN,
!C                                        SQRT
!C
!C***********************************************************
!C
      INTEGER LUN, NIT
      DOUBLE PRECISION A, A0, AA, B, B0, BND, C, C1, C2,                &
     &                 COSHM, COSHMM, D, D0, D1PD2, D2,                 &
     &                 DMAX, DSIG, DX, E, EMS, F, F0, FMAX,             &
     &                 FNEG, FTOL, R, RF, RSIG, RTOL, S, S1,            &
     &                 S2, SBIG, SCM, SIG, SINHM, SNEG,                 &
     &                 SSINH, SSM, STOL, T, T0, T1, T2, TM,             &
     &                 Y1L, Y2L
! in a module these need not to be declared here      DOUBLE PRECISION STORE
!C
      DATA SBIG/85.D0/,  LUN/-1/
!C
!C Store local parameters and test for errors.
!C
      RF = DBLE(IFL)
      DX = X2 - X1
      IF (DABS(RF) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 8
      Y1L = Y1
      Y2L = Y2
      BND = HBND
!C
!C Test for a valid constraint.
!C
      IF ((RF .LT. 0.D0  .AND.  MIN(Y1L,Y2L) .LT. BND)                  &
     &    .OR.  (RF .GT. 0.D0  .AND.                                    &
     &           BND .LT. MAX(Y1L,Y2L))) GO TO 9
!C
!C Test for infinite tension required.
!C
      S1 = Y1P
      S2 = Y2P
      IF ((Y1L .EQ. BND  .AND.  RF*S1 .GT. 0.D0)  .OR.                  &
     &    (Y2L .EQ. BND  .AND.  RF*S2 .LT. 0.D0)) GO TO 7
!C
!C Test for SIG = 0 sufficient.
!C
      SIG = 0.D0
      IF (RF*S1 .LE. 0.D0  .AND.  RF*S2 .GE. 0.D0) GO TO 6
!C
!C   Compute coefficients A0 and B0 of the Hermite cubic in-
!C     terpolant H0(x) = Y2 - DX*(S2*R + B0*R**2 + A0*R**3/3)
!C     where R = (X2-x)/DX.
!C
      S = (Y2L-Y1L)/DX
      T0 = 3.D0*S - S1 - S2
      A0 = 3.D0*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
!C
!C   H0 has local extrema in (X1,X2) iff S1*S2 < 0 or
!C     (T0*(S1+S2) < 0 and D0 .GE. 0).
!C
      IF (S1*S2 .GE. 0.D0  .AND.  (T0*(S1+S2) .GE. 0.D0                 &
     &    .OR.  D0 .LT. 0.D0)) GO TO 6
      IF (A0 .EQ. 0.D0) THEN
!C
!C   H0 is quadratic and has an extremum at R = -S2/(2*B0).
!C     H0(R) = Y2 + DX*S2**2/(4*B0).  Note that A0 = 0 im-
!C     plies 2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0.
!C     Also, the extremum is a min iff HBND is a lower bound.
!C
        F0 = (BND - Y2L - DX*S2*S2/(4.D0*B0))*RF
      ELSE
!C
!C   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
!C     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
!C     corresponds to a min.  The expression for R is chosen
!C     to avoid cancellation error.  H0(R) = Y2 + DX*(S2*B0 +
!C     2*D0*R)/(3*A0).
!C
        T = -B0 - SIGN(DSQRT(D0),B0)
        R = T/A0
        IF (RF*B0 .GT. 0.D0) R = S2/T
        F0 = (BND - Y2L - DX*(S2*B0+2.D0*D0*R)/(3.D0*A0))*RF
      ENDIF
!C
!C   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the
!C     constraint.
!C
      IF (F0 .GE. 0.D0) GO TO 6
!C
!C Find a zero of F(SIG) = (BND-H(R))*RF where the derivative
!C   of H, HP, vanishes at R.  F is a nondecreasing function,
!C   F(0) < 0, and F = FMAX for SIG sufficiently large.
!C
!C Initialize parameters for the secant method.  The method
!C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
!C   where SG0 and SNEG are defined implicitly by DSIG = SIG
!C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
!C   SNEG is initialized to a sufficiently large value that
!C   FNEG > 0.  This value is used only if the initial value
!C   of F is negative.
!C
      FMAX = MAX(1.D-3,MIN(dABS(Y1L-BND),dABS(Y2L-BND)))
      T = MAX(dABS(Y1L-BND),dABS(Y2L-BND))
      SIG = DX*MAX(dABS(S1),dABS(S2))/T
      DMAX = SIG*(1.D0-T/FMAX)
      SNEG = SIG - DMAX
      IF (LUN .GE. 0  .AND.  RF .LT. 0.D0)                              &
     &   WRITE (LUN,100) F0, FMAX, SNEG
      IF (LUN .GE. 0  .AND.  RF .GT. 0.D0)                              &
     &   WRITE (LUN,110) F0, FMAX, SNEG
  100 FORMAT (//1X,'SIG0 (LOWER BOUND) -- F(0) = ',D15.8,               &
     &        ', FMAX = ',D15.8/1X,46X,'SNEG = ',D15.8/)
  110 FORMAT (//1X,'SIG0 (UPPER BOUND) -- F(0) = ',D15.8,               &
     &        ', FMAX = ',D15.8/1X,46X,'SNEG = ',D15.8/)
      DSIG = SIG
      FNEG = FMAX
      D2 = S2 - S
      D1PD2 = S2 - S1
      NIT = 0
!C
!C Compute an absolute tolerance FTOL = abs(TOL), and a
!C   relative tolerance RTOL = 100*MACHEPS.
!C
      FTOL = dABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
!C
!C Top of loop:  compute F.
!C
    2 EMS = EXP(-SIG)
      IF (SIG .LE. .5D0) THEN
!C
!C   SIG .LE. .5:  use approximations designed to avoid can-
!C                 cellation error (associated with small
!C                 SIG) in the modified hyperbolic functions.
!C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        AA = A/EMS
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
!C
!C   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
!C              to avoid overflow.
!C
        TM = 1.D0 - EMS
        SSINH = TM*(1.D0+EMS)
        SSM = SSINH - 2.D0*SIG*EMS
        SCM = TM*TM
        C1 = SIG*SCM*D2 - SSM*D1PD2
        C2 = SIG*SSINH*D2 - SCM*D1PD2
        AA = 2.D0*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        A = EMS*AA
        E = SIG*SSINH - SCM - SCM
      ENDIF
!C
!C   HP(R) = S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E = 0
!C     for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D)),
!C     where ESR = EXP(SIG*R), A = C2-C1, D = B**2 - A*C, and
!C     B and C are defined below.
!C
      B = E*S2 - C2
      C = C2 + C1
      D = B*B - A*C
      F = 0.D0
      IF (AA*C .EQ. 0.D0  .AND.  B .EQ. 0.D0) GO TO 3
      F = FMAX
      IF (D .LT. 0.D0) GO TO 3
      T1 = DSQRT(D)
      T = -B - SIGN(T1,B)
      RSIG = 0.D0
      IF (RF*B .LT. 0.D0  .AND.  AA .NE. 0.) THEN
        IF (T/AA .GT. 0.D0) RSIG = SIG + LOG(T/AA)
      ENDIF
      IF ((RF*B .GT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.                   &
     &    C/T .GT. 0.D0) RSIG = LOG(C/T)
      IF ((RSIG .LE. 0.D0  .OR.  RSIG .GE. SIG)  .AND.                  &
     &    B .NE. 0.D0) GO TO 3
!C
!C   H(R) = Y2 - DX*(B*SIG*R + C1 + RF*SQRT(D))/(SIG*E).
!C
      F = (BND - Y2L + DX*(B*RSIG+C1+RF*T1)/(SIG*E))*RF
!C
!C   Update the number of iterations NIT.
!C
    3 NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120 FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',                     &
     &        D15.8)
      IF (F0*F .LT. 0.D0) THEN
!C
!C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
!C     and FNEG always have opposite signs.  If SIG is
!C     closer to SNEG than SG0, then swap (SNEG,FNEG) with
!C     (SG0,F0).
!C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF (dABS(DSIG) .GT. dABS(T1)) THEN
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
!C
!C   Test for convergence.
!C
      STOL = RTOL*SIG
      IF (dABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.               &
     &    F .LE. FTOL)  .OR.  dABS(F) .LE. RTOL) GO TO 6
!C
!C   Test for F0 = F = FMAX or F < 0 on the first iteration.
!C
      IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.  F .GT. 0.D0))            &
     &   GO TO 5
!C
!C   F*F0 > 0 and either the new estimate would be outside of
!C     the bracketing interval of length abs(DMAX) or F < 0
!C     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG).
!C
    4 DSIG = DMAX
      F0 = FNEG
!C
!C   Compute the change in SIG by linear interpolation be-
!C     tween (SG0,F0) and (SIG,F).
!C
    5 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130 FORMAT (1X,8X,'DSIG = ',D15.8)
      IF ( DABS(DSIG) .GT. DABS(DMAX)  .OR.                             &
     &     DSIG*DMAX .GT. 0. ) GO TO 4
!C
!C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
!C     Note that DSIG and DMAX have opposite signs.
!C
      IF (DABS(DSIG) .LT. STOL/2.D0)                                    &
     &  DSIG = -SIGN(STOL/2.D0,DMAX)
!C
!C   Bottom of loop:  update SIG, DMAX, and F0.
!C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 2
!C
!C No errors encountered and SIGMA finite.
!C
    6 IER = 0
      SIG0 = SIG
      RETURN
!C
!C Infinite tension required.
!C
    7 IER = 1
      SIG0 = SBIG
      RETURN
!C
!C Error in an input parameter.
!C
    8 IER = -1
      SIG0 = -1.D0
      RETURN
!C
!C Invalid constraint.
!C
    9 IER = -2
      SIG0 = -1.D0
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION SIG1 (X1,X2,Y1,Y2,Y1P,Y2P,              &
     &                                IFL,HPBND,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, HPBND, TOL
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/17/96
!C
!C   Given a pair of abscissae with associated ordinates and
!C slopes, this function determines the smallest (nonnega-
!C tive) tension factor SIGMA such that the derivative HP(x)
!C of the Hermite interpolatory tension spline H(x), defined
!C by SIGMA and the data, is bounded (either above or below)
!C by HPBND for all x in (X1,X2).
!C
!C On input:
!C
!C       X1,X2 = Abscissae.  X1 < X2.
!C
!C       Y1,Y2 = Values of H at X1 and X2.
!C
!C       Y1P,Y2P = Values of HP at X1 and X2.
!C
!C       IFL = Option indicator:
!C             IFL = -1 if HPBND is a lower bound on HP.
!C             IFL = 1 if HPBND is an upper bound on HP.
!C
!C       HPBND = Bound on HP.  If IFL = -1, HPBND .LE.
!C               min(Y1P,Y2P,S) for S = (Y2-Y1)/(X2-X1).  If
!C               IFL = 1, HPBND .GE. max(Y1P,Y2P,S).
!C
!C       TOL = Tolerance whose magnitude determines how close
!C             SIGMA is to its optimal value when nonzero
!C             finite tension is necessary and sufficient to
!C             satisfy the constraint.  For a lower bound,
!C             SIGMA is chosen so that HPBND .LE. HPMIN .LE.
!C             HPBND + abs(TOL), where HPMIN is the minimum
!C             value of HP in the interval, and for an upper
!C             bound, the maximum of HP satisfies HPBND -
!C             abs(TOL) .LE. HPMAX .LE. HPBND.  Thus, the
!C             constraint is satisfied but possibly with more
!C             tension than necessary.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered and the
!C                     constraint can be satisfied with fin-
!C                     ite tension.
!C             IER = 1 if no errors were encountered but in-
!C                     finite tension is required to satisfy
!C                     the constraint (e.g., IFL = -1, HPBND
!C                     = S, and Y1P > S).
!C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
!C             IER = -2 if HPBND is outside its valid range
!C                      on input.
!C
!C       SIG1 = Minimum tension factor defined above unless
!C              IER < 0, in which case SIG1 = -1.  If IER =
!C              1, SIG1 = 85, resulting in an approximation
!C              to the linear interpolant of the endpoint
!C              values.  Note, however, that SIG1 may be
!C              larger than 85 if IER = 0.
!C
!C Modules required by SIG1:  SNHCSH, STORE
!C
!C Intrinsic functions called by SIG1:  ABS, DBLE, EXP, MAX,
!C                                        MIN, SIGN, SQRT
!C
!C***********************************************************
!C
      INTEGER LUN, NIT
      DOUBLE PRECISION A, A0, B0, BND, C0, C1, C2, COSHM,               &
     &                 COSHMM, D0, D1, D1PD2, D2, DMAX,                 &
     &                 DSIG, DX, E, EMS, EMS2, F, F0, FMAX,             &
     &                 FNEG, FTOL, RF, RTOL, S, S1, S2,                 &
     &                 SBIG, SIG, SINH, SINHM, STOL, T0, T1,            &
     &                 T2, TM
! in a module these need not to be declared here      DOUBLE PRECISION STORE
!C
      DATA SBIG/85.D0/,  LUN/-1/
!C
!C Store local parameters and test for errors.
!C
      RF = DBLE(IFL)
      DX = X2 - X1
      IF (DABS(RF) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 7
      S1 = Y1P
      S2 = Y2P
      S = (Y2-Y1)/DX
      BND = HPBND
!C
!C Test for a valid constraint.
!C
      IF ((RF .LT. 0.D0  .AND.  MIN(S1,S2,S) .LT. BND)                  &
     &    .OR.  (RF .GT. 0.D0  .AND.                                    &
     &           BND .LT. MAX(S1,S2,S))) GO TO 8
!C
!C Test for infinite tension required.
!C
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))              &
     &   GO TO 6
!C
!C Test for SIG = 0 sufficient.  The Hermite cubic interpo-
!C   land H0 has derivative HP0(x) = S2 + 2*B0*R + A0*R**2,
!C   where R = (X2-x)/DX.
!C
      SIG = 0.D0
      T0 = 3.D0*S - S1 - S2
      B0 = T0 - S2
      C0 = T0 - S1
      A0 = -B0 - C0
!C
!C   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
!C     B0*C0 > 0 and the third derivative of H0 has the
!C     sign of A0.
!C
      IF (B0*C0 .LE. 0.D0  .OR.  A0*RF .GT. 0.D0) GO TO 5
!C
!C   A0*RF < 0 and HP0(R) = -D0/A0 at R = -B0/A0.
!C
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/A0)*RF
      IF (F0 .GE. 0.D0) GO TO 5
!C
!C Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an
!C   extremum at R.  F has a unique zero, F(0) = F0 < 0, and
!C   F = (BND-S)*RF > 0 for SIG sufficiently large.
!C
!C Initialize parameters for the secant method.  The method
!C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
!C   where SG0 and SNEG are defined implicitly by DSIG = SIG
!C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
!C   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/
!C   (SIG-2.))*RF -- a value for which F(SIG) .GE. 0 and
!C   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
!C   significant relative to EXP(SIG).
!C
      FMAX = (BND-S)*RF
      IF (LUN .GE. 0  .AND.  RF .LT. 0.D0)                              &
     &  WRITE (LUN,100) F0, FMAX
      IF (LUN .GE. 0  .AND.  RF .GT. 0.D0)                              &
     &  WRITE (LUN,110) F0, FMAX
  100 FORMAT (//1X,'SIG1 (LOWER BOUND) -- F(0) = ',D15.8,               &
     &        ', FMAX = ',D15.8/)
  110 FORMAT (//1X,'SIG1 (UPPER BOUND) -- F(0) = ',D15.8,               &
     &        ', FMAX = ',D15.8/)
      SIG = 2.D0 - A0/(3.D0*(BND-S))
      IF (STORE(SIG*EXP(-SIG)+.5D0) .EQ. .5D0) GO TO 5
      DSIG = SIG
      DMAX = -2.D0*SIG
      FNEG = FMAX
      D1 = S - S1
      D2 = S2 - S
      D1PD2 = D1 + D2
      NIT = 0
!C
!C Compute an absolute tolerance FTOL = abs(TOL), and a
!C   relative tolerance RTOL = 100*MACHEPS.
!C
      FTOL = DABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
!C
!C Top of loop:  compute F.
!C
    2 IF (SIG .LE. .5D0) THEN
!C
!C   Use approximations designed to avoid cancellation error
!C     (associated with small SIG) in the modified hyperbolic
!C     functions.
!C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
!C
!C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
!C     overflow.
!C
        EMS = EXP(-SIG)
        EMS2 = EMS + EMS
        TM = 1.D0 - EMS
        SINH = TM*(1.D0+EMS)
        SINHM = SINH - SIG*EMS2
        COSHM = TM*TM
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*SINH*D2 - COSHM*D1PD2
        A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        E = SIG*SINH - COSHM - COSHM
      ENDIF
!C
!C   The second derivative of H(R) has a zero at EXP(SIG*R) =
!C     SQRT((C2+C1)/A) and R is in (0,1) and well-defined
!C     iff HPP(X1)*HPP(X2) < 0.
!C
      F = FMAX
      T1 = A*(C2+C1)
      IF (T1 .GE. 0.D0) THEN
        IF (C1*(SIG*COSHM*D1 - SINHM*D1PD2) .LT. 0.D0) THEN
!C
!C   HP(R) = (B+SIGN(A)*SQRT(A*!C))/E at the critical value
!C     of R, where A = C2-C1, B = E*S2-C2, and !C = C2+C1.
!C     NOTE THAT RF*A < 0.
!C
          F = (BND - (E*S2-C2 - RF*DSQRT(T1))/E)*RF
        ENDIF
      ENDIF
!C
!C   Update the number of iterations NIT.
!C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120 FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',                     &
     &        D15.8)
      IF (F0*F .LT. 0.D0) THEN
!C
!C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
!C     and FNEG always have opposite signs.  If SIG is closer
!C     to SNEG than SG0 and abs(F) < abs(FNEG), then swap
!C     (SNEG,FNEG) with (SG0,F0).
!C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF ( DABS(DSIG) .GT. DABS(T1)  .AND.                            &
     &       DABS(F) .LT. DABS(T2) ) THEN
!C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
!C
!C   Test for convergence.
!C
      STOL = RTOL*SIG
      IF (DABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.               &
     &    F .LE. FTOL)  .OR.  DABS(F) .LE. RTOL) GO TO 5
      IF (F0*F .LT. 0  .OR.  DABS(F) .LT. DABS(F0)) GO TO 4
!C
!C   F*F0 > 0 and the new estimate would be outside of the
!C     bracketing interval of length abs(DMAX).  Reset
!C     (SG0,F0) to (SNEG,FNEG).
!C
    3 DSIG = DMAX
      F0 = FNEG
!C
!C   Compute the change in SIG by linear interpolation be-
!C     tween (SG0,F0) and (SIG,F).
!C
    4 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130 FORMAT (1X,8X,'DSIG = ',D15.8)
      IF ( DABS(DSIG) .GT. DABS(DMAX)  .OR.                             &
     &     DSIG*DMAX .GT. 0. ) GO TO 3
!C
!C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
!C     Note that DSIG and DMAX have opposite signs.
!C
      IF (DABS(DSIG) .LT. STOL/2.D0)                                    &
     &  DSIG = -SIGN(STOL/2.D0,DMAX)
!C
!C   Bottom of loop:  update SIG, DMAX, and F0.
!C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 2
!C
!C No errors encountered and SIGMA finite.
!C
    5 IER = 0
      SIG1 = SIG
      RETURN
!C
!C Infinite tension required.
!C
    6 IER = 1
      SIG1 = SBIG
      RETURN
!C
!C Error in an input parameter.
!C
    7 IER = -1
      SIG1 = -1.D0
      RETURN
!C
!C Invalid constraint.
!C
    8 IER = -2
      SIG1 = -1.D0
      RETURN
      END FUNCTION
      
      DOUBLE PRECISION FUNCTION SIG2 (X1,X2,Y1,Y2,Y1P,Y2P,              &
     &                                IFL,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, TOL
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   07/08/92
!C
!C   Given a pair of abscissae with associated ordinates and
!C slopes, this function determines the smallest (nonnega-
!C tive) tension factor SIGMA such that the Hermite interpo-
!C latory tension spline H(x) preserves convexity (or con-
!C cavity) of the data;  i.e.,
!C
!C   Y1P .LE. S .LE. Y2P implies HPP(x) .GE. 0  or
!C   Y1P .GE. S .GE. Y2P implies HPP(x) .LE. 0
!C
!C for all x in the open interval (X1,X2), where S = (Y2-Y1)/
!C (X2-X1) and HPP denotes the second derivative of H.  Note,
!C however, that infinite tension is required if Y1P = S or
!C Y2P = S (unless Y1P = Y2P = S).
!C
!C On input:
!C
!C       X1,X2 = Abscissae.  X1 < X2.
!C
!C       Y1,Y2 = Values of H at X1 and X2.
!C
!C       Y1P,Y2P = Derivative values of H at X1 and X2.
!C
!C       IFL = Option indicator (sign of HPP):
!C             IFL = -1 if HPP is to be bounded above by 0.
!C             IFL = 1 if HPP is to be bounded below by 0
!C                     (preserve convexity of the data).
!C
!C       TOL = Tolerance whose magnitude determines how close
!C             SIGMA is to its optimal value when nonzero
!C             finite tension is necessary and sufficient to
!C             satisfy convexity or concavity.  In the case
!C             of convexity, SIGMA is chosen so that 0 .LE.
!C             HPPMIN .LE. abs(TOL), where HPPMIN is the min-
!C             imum value of HPP in the interval.  In the
!C             case of concavity, the maximum value of HPP
!C             satisfies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus,
!C             the constraint is satisfied but possibly with
!C             more tension than necessary.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered and fin-
!C                     ite tension is sufficient to satisfy
!C                     the constraint.
!C             IER = 1 if no errors were encountered but in-
!C                     finite tension is required to satisfy
!C                     the constraint.
!C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
!C             IER = -2 if the constraint cannot be satis-
!C                      fied:  the sign of S-Y1P or Y2P-S
!C                      does not agree with IFL.
!C
!C       SIG2 = Tension factor defined above unless IER < 0,
!C              in which case SIG2 = -1.  If IER = 1, SIG2
!C              is set to 85, resulting in an approximation
!C              to the linear interpolant of the endpoint
!C              values.  Note, however, that SIG2 may be
!C              larger than 85 if IER = 0.
!C
!C Modules required by SIG2:  SNHCSH, STORE
!C
!C Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
!C                                        SQRT
!C
!C***********************************************************
!C
      INTEGER LUN, NIT
      DOUBLE PRECISION COSHM, D1, D2, DSIG, DUMMY, DX, EMS,             &
     &                 F, FP, FTOL, RTOL, S, SBIG, SIG,                 &
     &                 SINHM, SSM, T, T1, TP1
! in a module these need not to be declared here      DOUBLE PRECISION STORE
!C
      DATA SBIG/85.D0/, LUN/-1/
!C
!C Test for an errors in the input parameters.
!C
      DX = X2 - X1
!      IF (ABS(IFL) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 5
      IF (ABS(IFL) .NE. 1 .OR.  DX .LE. 0.D0) GO TO 5
!C
!C Compute the slope and second differences, and test for
!C   an invalid constraint.
!C
      S = (Y2-Y1)/DX
      D1 = S - Y1P
      D2 = Y2P - S
      IF ((IFL .GT. 0.D0  .AND.  MIN(D1,D2) .LT. 0.D0)                  &
     &    .OR.  (IFL .LT. 0.D0  .AND.                                   &
     &           MAX(D1,D2) .GT. 0.D0)) GO TO 6
!C
!C Test for infinite tension required.
!C
      IF (D1*D2 .EQ. 0.D0  .AND.  D1 .NE. D2) GO TO 4
!C
!C Test for SIG = 0 sufficient.
!C
      SIG = 0.D0
      IF (D1*D2 .EQ. 0.D0) GO TO 3
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.D0) GO TO 3
!C
!C Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
!C   Since the derivative of F vanishes at the origin, a
!C   quadratic approximation is used to obtain an initial
!C   estimate for the Newton method.
!C
      TP1 = T + 1.D0
      SIG = DSQRT(10.D0*T-20.D0)
      NIT = 0
!C
!C   Compute an absolute tolerance FTOL = abs(TOL) and a
!C     relative tolerance RTOL = 100*MACHEPS.
!C
      FTOL = DABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
!C
!C Evaluate F and its derivative FP.
!C
    2 IF (SIG .LE. .5D0) THEN
!C
!C   Use approximations designed to avoid cancellation error
!C     in the hyperbolic functions.
!C
        CALL SNHCSH (SIG, SINHM,COSHM,DUMMY)
        T1 = COSHM/SINHM
        FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.D0)
      ELSE
!C
!C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
!C     overflow.
!C
        EMS = EXP(-SIG)
        SSM = 1.D0 - EMS*(EMS+SIG+SIG)
        T1 = (1.D0-EMS)*(1.D0-EMS)/SSM
        FP = T1 + SIG*(2.D0*SIG*EMS/SSM - T1*T1 + 1.D0)
      ENDIF
!C
      F = SIG*T1 - TP1
      IF (LUN .GE. 0) WRITE (LUN,100) SIG, F, FP
  100 FORMAT (1X,'SIG2 -- SIG = ',D15.8,', F(SIG) = ',                  &
     &        D15.8/1X,29X,'FP(SIG) = ',D15.8)
      NIT = NIT + 1
!C
!C   Test for convergence.
!C
      IF (FP .LE. 0.D0) GO TO 3
      DSIG = -F/FP
      IF (DABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.D0 .AND.            &
     &    F .LE. FTOL)  .OR.  DABS(F) .LE. RTOL) GO TO 3
!C
!C   Update SIG.
!C
      SIG = SIG + DSIG
      GO TO 2
!C
!C No errors encountered, and SIGMA is finite.
!C
    3 IER = 0
      SIG2 = SIG
      RETURN
!C
!C Infinite tension required.
!C
    4 IER = 1
      SIG2 = SBIG
      RETURN
!C
!C X2 .LE. X1 or abs(IFL) .NE. 1.
!C
    5 IER = -1
      SIG2 = -1.D0
      RETURN
!C
!C The constraint cannot be satisfied.
!C
    6 IER = -2
      SIG2 = -1.D0
      RETURN
      END FUNCTION
      
      SUBROUTINE SIGBI (N,X,Y,YP,TOL,B,BMAX, SIGMA, ICFLG,              &
     &                  DSMAX,IER)
      INTEGER N, ICFLG(N), IER
      DOUBLE PRECISION X(N), Y(N), YP(N), TOL, B(5,N), BMAX,            &
     &                 SIGMA(N), DSMAX
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   Given a set of abscissae X with associated data values Y
!C and derivatives YP, this subroutine determines the small-
!C est (nonnegative) tension factors SIGMA such that the Her-
!C mite interpolatory tension spline H(x) satisfies a set of
!C user-specified constraints.
!C
!C   SIGBI may be used in conjunction with Subroutine YPC2
!C (or YPC2P) in order to produce a !C-2 interpolant which
!C satisfies the constraints.  This is achieved by calling
!C YPC2 with SIGMA initialized to the zero vector, and then
!C alternating calls to SIGBI with calls to YPC2 until the
!C change in SIGMA is small (refer to the parameter descrip-
!C tions for SIGMA, DSMAX and IER), or the maximum relative
!C change in YP is bounded by a tolerance (a reasonable value
!C is .01).  A similar procedure may be used to produce a !C-2
!C shape-preserving smoothing curve (Subroutine SMCRV).
!C
!C   Refer to Subroutine SIGS for a means of selecting mini-
!C mum tension factors to preserve shape properties of the
!C data.
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing a strictly in-
!C           creasing sequence of abscissae:  X(I) < X(I+1)
!C           for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values (or
!C           function values computed by SMCRV) associated
!C           with the abscissae.  H(X(I)) = Y(I) for I =
!C           1,...,N.
!C
!C       YP = Array of length N containing first derivatives
!C            of H at the abscissae.  Refer to Subroutines
!C            YPC1, YPC1P, YPC2, YPC2P, and SMCRV.
!C
!C       TOL = Tolerance whose magnitude determines how close
!C             each tension factor is to its optimal value
!C             when nonzero finite tension is necessary and
!C             sufficient to satisfy a constraint.  Refer to
!C             functions SIG0, SIG1, and SIG2.  TOL should be
!C             set to 0 for optimal tension.
!C
!C       B = Array dimensioned 5 by N-1 containing bounds or
!C           flags which define the constraints.  For I = 1
!C           to N-1, column I defines the constraints associ-
!C           ated with interval I (X(I),X(I+1)) as follows:
!C
!C             B(1,I) is an upper bound on H
!C             B(2,I) is a lower bound on H
!C             B(3,I) is an upper bound on HP
!C             B(4,I) is a lower bound on HP
!C             B(5,I) specifies the required sign of HPP
!C
!C           where HP and HPP denote the first and second
!C           derivatives of H, respectively.  A null con-
!C           straint is specified by abs(B(K,I)) .GE. BMAX
!C           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX,
!C           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE.
!C           -BMAX, or B(5,I) = 0.  Any positive value of
!C           B(5,I) specifies that H should be convex, a
!C           negative values specifies that H should be con-
!C           cave, and 0 specifies that no restriction be
!C           placed on HPP.  Refer to Functions SIG0, SIG1,
!C           and SIG2 for definitions of valid constraints.
!C
!C       BMAX = User-defined value of infinity which, when
!C              used as an upper bound in B (or when its
!C              negative is used as a lower bound), specifies
!C              that no constraint is to be enforced.
!C
!C The above parameters are not altered by this routine.
!C
!C       SIGMA = Array of length N-1 containing minimum val-
!C               ues of the tension factors.  SIGMA(I) is as-
!C               sociated with interval (I,I+1) and SIGMA(I)
!C               .GE. 0 for I = 1,...,N-1.  SIGMA should be
!C               set to the zero vector if minimal tension
!C               is desired, and should be unchanged from a
!C               previous call in order to ensure convergence
!C               of the !C-2 iterative procedure.
!C
!C       ICFLG = Array of length .GE. N-1.
!C
!C On output:
!C
!C       SIGMA = Array containing tension factors for which
!C               H(x) satisfies the constraints defined by B,
!C               with the restriction that SIGMA(I) .LE. 85
!C               for all I (unless the input value is larger).
!C               The factors are as small as possible (within
!C               the tolerance), but not less than their
!C               input values.  If infinite tension is re-
!C               quired in interval (X(I),X(I+1)), then
!C               SIGMA(I) = 85 (and H is an approximation to
!C               the linear interpolant on the interval),
!C               and if no constraint is specified in the
!C               interval, then SIGMA(I) = 0 (unless the
!C               input value is positive), and thus H is
!C               cubic.  Invalid constraints are treated as
!C               null constraints.
!C
!C       ICFLG = Array of invalid constraint flags associated
!C               with intervals.  For I = 1 to N-1, ICFLG(I)
!C               is a 5-bit value b5b4b3b2b1, where bK = 1 if
!C               and only if constraint K cannot be satis-
!C               fied.  Thus, all constraints in interval I
!C               are satisfied if and only if ICFLG(I) = 0
!C               (and IER .GE. 0).
!C
!C       DSMAX = Maximum increase in a component of SIGMA
!C               from its input value.  The increase is a
!C               relative change if the input value is
!C               positive, and an absolute change otherwise.
!C
!C       IER = Error indicator and information flag:
!C             IER = I if no errors (other than invalid con-
!C                     straints) were encountered and I
!C                     components of SIGMA were altered from
!C                     their input values for 0 .LE. I .LE.
!C                     N-1.
!C             IER = -1 if N < 2.  SIGMA and ICFLG are not
!C                      altered in this case.
!C             IER = -I if X(I) .LE. X(I-1) for some I in the
!C                      range 2,...,N.  SIGMA(J) and ICFLG(J)
!C                      are unaltered for J .GE. I-1 in this
!C                      case.
!C
!C Modules required by SIGBI:  SIG0, SIG1, SIG2, SNHCSH,
!C                               STORE
!C
!C Intrinsic functions called by SIGBI:  ABS, MAX, MIN
!C
!C***********************************************************
!C
      INTEGER I, ICFK, ICNT, IERR, IFL, K, NM1
      DOUBLE PRECISION BMX, BND, DSIG, DSM, S, SBIG, SIG,               &
     &                 SIGIN
! in a module these need not to be declared here      DOUBLE PRECISION SIG0, SIG1, SIG2
!C
      DATA SBIG/85.D0/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 4
      BMX = BMAX
!C
!C Initialize change counter ICNT and maximum change DSM for
!C   loop on intervals.
!C
      ICNT = 0
      DSM = 0.D0
      DO 3 I = 1,NM1
        IF (X(I) .GE. X(I+1)) GO TO 5
        ICFLG(I) = 0
!C
!C Loop on constraints for interval I.  SIG is set to the
!C   largest tension factor required to satisfy all five
!C   constraints.  ICFK = 2**(K-1) is the increment for
!C   ICFLG(I) when constraint K is invalid.
!C
        SIG = 0.D0
        ICFK = 1
        DO 2 K = 1,5
          BND = B(K,I)
          IF (K .LT. 5  .AND.  DABS(BND) .GE. BMX) GO TO 1
          IF (K .LE. 2) THEN
            IFL = 3 - 2*K
            S = SIG0 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),            &
     &                IFL,BND,TOL, IERR)
          ELSEIF (K .LE. 4) THEN
            IFL = 7 - 2*K
            S = SIG1 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),            &
     &                IFL,BND,TOL, IERR)
          ELSE
            IF (BND .EQ. 0.D0) GO TO 1
            IFL = -1
            IF (BND .GT. 0.D0) IFL = 1
            S = SIG2 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),            &
     &                IFL,TOL, IERR)
          ENDIF
          IF (IERR .EQ. -2) THEN
!C
!C   An invalid constraint was encountered.  Increment
!C     ICFLG(I).
!C
            ICFLG(I) = ICFLG(I) + ICFK
          ELSE
!C
!C   Update SIG.
!C
            SIG = MAX(SIG,S)
          ENDIF
!C
!C   Bottom of loop on constraints K:  update ICFK.
!C
    1     ICFK = 2*ICFK
    2     CONTINUE
!C
!C Bottom of loop on intervals:  update SIGMA(I), ICNT, and
!C   DSM if necessary.
!C
        SIG = MIN(SIG,SBIG)
        SIGIN = SIGMA(I)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    3   CONTINUE
!C
!C No errors (other than invalid constraints) encountered.
!C
      DSMAX = DSM
      IER = ICNT
      RETURN
!C
!C N < 2.
!C
    4 DSMAX = 0.D0
      IER = -1
      RETURN
!C
!C X(I) .GE. X(I+1).
!C
    5 DSMAX = DSM
      IER = -(I+1)
      RETURN
      END SUBROUTINE
      
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      DOUBLE PRECISION X, SINHM, COSHM, COSHMM
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/20/96
!C
!C   This subroutine computes approximations to the modified
!C hyperbolic functions defined below with relative error
!C bounded by 3.4E-20 for a floating point number system with
!C sufficient precision.
!C
!C   Note that the 21-digit constants in the data statements
!C below may not be acceptable to all compilers.
!C
!C On input:
!C
!C       X = Point at which the functions are to be
!C           evaluated.
!C
!C X is not altered by this routine.
!C
!C On output:
!C
!C       SINHM = sinh(X) - X.
!C
!C       COSHM = cosh(X) - 1.
!C
!C       COSHMM = cosh(X) - 1 - X*X/2.
!C
!C Modules required by SNHCSH:  None
!C
!C Intrinsic functions called by SNHCSH:  ABS, EXP
!C
!C***********************************************************
!C
      DOUBLE PRECISION AX, EXPX, F, P, P1, P2, P3, P4, Q,               &
     &                 Q1, Q2, Q3, Q4, XC, XS, XSD2, XSD4
!C
      DATA P1/-3.51754964808151394800D5/,                               &
     &     P2/-1.15614435765005216044D4/,                               &
     &     P3/-1.63725857525983828727D2/,                               &
     &     P4/-7.89474443963537015605D-1/
      DATA Q1/-2.11052978884890840399D6/,                               &
     &     Q2/3.61578279834431989373D4/,                                &
     &     Q3/-2.77711081420602794433D2/,                               &
     &     Q4/1.D0/
      AX = DABS(X)
      XS = AX*AX
      IF (AX .LE. .5D0) THEN
!C
!C Approximations for small X:
!C
        XC = X*XS
        P = ((P4*XS+P3)*XS+P2)*XS+P1
        Q = ((Q4*XS+Q3)*XS+Q2)*XS+Q1
        SINHM = XC*(P/Q)
        XSD4 = .25D0*XS
        XSD2 = XSD4 + XSD4
        P = ((P4*XSD4+P3)*XSD4+P2)*XSD4+P1
        Q = ((Q4*XSD4+Q3)*XSD4+Q2)*XSD4+Q1
        F = XSD4*(P/Q)
        COSHMM = XSD2*F*(F+2.D0)
        COSHM = COSHMM + XSD2
      ELSE
!C
!C Approximations for large X:
!C
        EXPX = EXP(AX)
        SINHM = -(((1.D0/EXPX+AX)+AX)-EXPX)/2.D0
        IF (X .LT. 0.D0) SINHM = -SINHM
        COSHM = ((1.D0/EXPX-2.D0)+EXPX)/2.D0
        COSHMM = COSHM - XS/2.D0
      ENDIF
      RETURN
      END SUBROUTINE
      
      DOUBLE PRECISION FUNCTION STORE (X)
      DOUBLE PRECISION X
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This function forces its argument X to be stored in a
!C memory location, thus providing a means of determining
!C floating point number characteristics (such as the machine
!C precision) when it is necessary to avoid computation in
!C high precision registers.
!C
!C On input:
!C
!C       X = Value to be stored.
!C
!C X is not altered by this function.
!C
!C On output:
!C
!C       STORE = Value of X after it has been stored and
!C               possibly truncated or rounded to the single
!C               precision word length.
!C
!C Modules required by STORE:  None
!C
!C***********************************************************
!C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
      Y = X
      STORE = Y
      RETURN
      END FUNCTION
      
      SUBROUTINE TSPBI (N,X,Y,NCD,IENDC,PER,B,BMAX,LWK, WK,             &
     &                  YP, SIGMA,ICFLG,IER)
      INTEGER N, NCD, IENDC, LWK, ICFLG(N), IER
      LOGICAL PER
      DOUBLE PRECISION X(N), Y(N), B(5,N), BMAX, WK(LWK),               &
     &                 YP(N), SIGMA(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   07/08/92
!C
!C   This subroutine computes a set of parameter values which
!C define a Hermite interpolatory tension spline H(x).  The
!C parameters consist of knot derivative values YP computed
!C by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension
!C factors SIGMA chosen to satisfy user-specified constraints
!C (by Subroutine SIGBI).  Refer to Subroutine TSPSI for an
!C alternative method of computing tension factors.
!C
!C   Refer to Subroutine TSPSS for a means of computing
!C parameters which define a smoothing curve rather than an
!C interpolatory curve.
!C
!C   The tension spline may be evaluated by Subroutine TSVAL1
!C or Functions HVAL (values), HPVAL (first derivatives),
!C HPPVAL (second derivatives), and TSINTL (integrals).
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 2 and N .GE. 3 if
!C           PER = TRUE.
!C
!C       X = Array of length N containing a strictly in-
!C           creasing sequence of abscissae:  X(I) < X(I+1)
!C           for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values asso-
!C           ciated with the abscissae.  H(X(I)) = Y(I) for
!C           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N)
!C           is set to Y(1).
!C
!C       NCD = Number of continuous derivatives at the knots.
!C             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
!C             are computed by local monotonicity-constrained
!C             quadratic fits.  Otherwise, a linear system is
!C             solved for the derivative values which result
!C             in second derivative continuity.  This re-
!C             quires iterating on calls to YPC2 or YPC2P and
!C             calls to SIGBI, and generally results in more
!C             nonzero tension factors (hence more expensive
!C             evaluation).
!C
!C       IENDC = End condition indicator for NCD = 2 and PER
!C               = FALSE (or dummy parameter otherwise):
!C               IENDC = 0 if YP(1) and YP(N) are to be com-
!C                         puted by monotonicity-constrained
!C                         parabolic fits to the first three
!C                         and last three points, respective-
!C                         ly.  This is identical to the
!C                         values computed by YPC1.
!C               IENDC = 1 if the first derivatives of H at
!C                         X(1) and X(N) are user-specified
!C                         in YP(1) and YP(N), respectively.
!C               IENDC = 2 if the second derivatives of H at
!C                         X(1) and X(N) are user-specified
!C                         in YP(1) and YP(N), respectively.
!C               IENDC = 3 if the end conditions are to be
!C                         computed by Subroutine ENDSLP and
!C                         vary with SIGMA(1) and SIGMA(N-1).
!C
!C       PER = Logical variable with value TRUE if and only
!C             H(x) is to be a periodic function with period
!C             X(N)-X(1).  It is assumed without a test that
!C             Y(N) = Y(1) in this case.  On output, YP(N) =
!C             YP(1).  If H(x) is one of the components of a
!C             parametric curve, this option may be used to
!C             obtained a closed curve.
!C
!C       B = Array dimensioned 5 by N-1 containing bounds or
!C           flags which define the constraints.  For I = 1
!C           to N-1, column I defines the constraints associ-
!C           ated with interval (X(I),X(I+1)) as follows:
!C
!C             B(1,I) is an upper bound on H
!C             B(2,I) is a lower bound on H
!C             B(3,I) is an upper bound on HP
!C             B(4,I) is a lower bound on HP
!C             B(5,I) specifies the required sign of HPP
!C
!C           where HP and HPP denote the first and second
!C           derivatives of H, respectively.  A null con-
!C           straint is specified by abs(B(K,I)) .GE. BMAX
!C           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX,
!C           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE.
!C           -BMAX, or B(5,I) = 0.  Any positive value of
!C           B(5,I) specifies that H should be convex, a
!C           negative values specifies that H should be con-
!C           cave, and 0 specifies that no restriction be
!C           placed on HPP.  Refer to Functions SIG0, SIG1,
!C           and SIG2 for definitions of valid constraints.
!C
!C       BMAX = User-defined value of infinity which, when
!C              used as an upper bound in B (or when when
!C              its negative is used as a lower bound),
!C              specifies that no constraint is to be en-
!C              forced.
!C
!C       LWK = Length of work space WK:
!C             LWK GE 2N-2 if NCD = 2 and PER = FALSE
!C             LWK GE 3N-3 if NCD = 2 and PER = TRUE
!C
!C   The above parameters, except possibly Y(N), are not
!C altered by this routine.
!C
!C       WK = Array of length at least LWK to be used as
!C            temporary work space.
!C
!C       YP = Array of length .GE. N containing end condition
!C            values in positions 1 and N if NCD = 2 and
!C            IENDC = 1 or IENDC = 2.
!C
!C       SIGMA = Array of length .GE. N-1.
!C
!C       ICFLG = Array of length .GE. N-1.
!C
!C On output:
!C
!C       WK = Array containing convergence parameters in the
!C            first two locations if IER > 0 (NCD = 2 and
!C            no error other than invalid constraints was
!C            encountered):
!C            WK(1) = Maximum relative change in a component
!C                    of YP on the last iteration.
!C            WK(2) = Maximum relative change in a component
!C                    of SIGMA on the last iteration.
!C
!C       YP = Array containing derivatives of H at the
!C            abscissae.  YP is not altered if -3 < IER < 0,
!C            and YP is only partially defined if IER = -4.
!C
!C       SIGMA = Array containing tension factors for which
!C               H(x) satisfies the constraints defined by B.
!C               SIGMA(I) is associated with interval (X(I),
!C               X(I+1)) for I = 1,...,N-1.  If infinite ten-
!C               sion is required in interval I, then
!C               SIGMA(I) = 85 (and H is an approximation to
!C               the linear interpolant on the interval),
!C               and if no constraint is specified in the
!C               interval, then SIGMA(I) = 0, and thus H is
!C               cubic.  Invalid constraints are treated as
!C               null constraints.  SIGMA is not altered if
!C               -3 < IER < 0 (unless IENDC is invalid), and
!C               SIGMA is the zero vector if IER = -4 or
!C               IENDC (if used) is invalid.
!C
!C       ICFLG = Array of invalid constraint flags associated
!C               with intervals.  For I = 1 to N-1, ICFLG(I)
!C               is a 5-bit value b5b4b3b2b1, where bK = 1 if
!C               and only if constraint K cannot be satis-
!C               fied.  Thus, all constraints in interval I
!C               are satisfied if and only if ICFLG(I) = 0
!C               (and IER .GE. 0).  ICFLG is not altered if
!C               IER < 0.
!C
!C       IER = Error indicator or iteration count:
!C             IER = IC .GE. 0 if no errors were encountered
!C                      (other than invalid constraints) and
!C                      IC calls to SIGBI and IC+1 calls to
!C                      YPC1, YPC1P, YPC2 or YPC2P were
!C                      employed.  (IC = 0 if NCD = 1).
!C             IER = -1 if N, NCD, or IENDC is outside its
!C                      valid range.
!C             IER = -2 if LWK is too small.
!C             IER = -4 if the abscissae X are not strictly
!C                      increasing.
!C
!C Modules required by TSPBI:  ENDSLP, SIG0, SIG1, SIG2,
!C                               SIGBI, SNHCSH, STORE,
!C                               YPCOEF, YPC1, YPC1P, YPC2,
!C                               YPC2P
!C
!C Intrinsic functions called by TSPBI:  ABS, MAX
!C
!C***********************************************************
!C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      LOGICAL LOOP2
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, E, STOL, YP1, YPN
!C
      DATA STOL/0.D0/,  MAXIT/49/,  DYPTOL/.01D0/
!C
!C Convergence parameters:
!C
!C   STOL = Absolute tolerance for SIGBI.
!C   MAXIT = Maximum number of YPC2/SIGBI iterations for
!C             each loop if NCD = 2.
!C   DYPTOL = Bound on the maximum relative change in a
!C              component of YP defining convergence of
!C              the YPC2/SIGBI iteration when NCD = 2.
!C
      NN = N
      NM1 = NN - 1
!C
!C Test for invalid input parameters N, NCD, or LWK.
!C
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.                &
     &    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF ( NCD .EQ. 2  .AND.  (LWK .LT. 2*NM1  .OR.                     &
     &     (PER  .AND.  LWK .LT. 3*NM1)) ) GO TO 12
!C
!C Initialize iteration count ITER, and initialize SIGMA to
!C   zeros.
!C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = 0.D0
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
!C
!C NCD = 1.
!C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,X,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,X,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 14
!C
!C   Compute tension factors.
!C
        CALL SIGBI (NN,X,Y,YP,STOL,B,BMAX, SIGMA, ICFLG,                &
     &              DSMAX,IERR)
        GO TO 10
      ENDIF
!C
!C NCD = 2.
!C
      IF (.NOT. PER) THEN
!C
!C   Nonperiodic case:  call YPC2 and test for IENDC or X
!C     invalid.
!C
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,                    &
     &             WK, YP,IERR)
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 14
      ELSE
!C
!C   Periodic fit:  call YPC2P.
!C
        CALL YPC2P (NN,X,Y,SIGMA,WK, YP,IERR)
        IF (IERR .GT. 1) GO TO 14
      ENDIF
      LOOP2 = .FALSE.
!C
!C   Iterate on calls to SIGBI and YPC2 (or YPC2P).  The
!C     first N-1 WK locations are used to store the deriva-
!C     tive estimates YP from the previous iteration.
!C
!C   LOOP2 is TRUE iff tension factors are not allowed to
!C         decrease between iterations (loop 1 failed to
!C         converge with MAXIT iterations).
!C   DYP is the maximum relative change in a component of YP.
!C   ICNT is the number of tension factors which were altered
!C        by SIGBI.
!C   DSMAX is the maximum relative change in a component of
!C         SIGMA.
!C
    2 DO 6 ITER = 1,MAXIT
        DYP = 0.D0
        DO 3 I = 2,NM1
          WK(I) = YP(I)
    3     CONTINUE
        CALL SIGBI (NN,X,Y,YP,STOL,B,BMAX, SIGMA, ICFLG,                &
     &              DSMAX,ICNT)
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,                  &
     &               WK(NN), YP,IERR)
        ELSE
          CALL YPC2P (NN,X,Y,SIGMA,WK(NN), YP,IERR)
        ENDIF
        DO 4 I = 2,NM1
          E = DABS(YP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/DABS(WK(I))
          DYP = MAX(DYP,E)
    4     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 7
        IF (.NOT. LOOP2) THEN
!C
!C   Loop 1:  reinitialize SIGMA to zeros.
!C
          DO 5 I = 1,NM1
            SIGMA(I) = 0.D0
    5       CONTINUE
        ENDIF
    6   CONTINUE
!C
!C   The loop failed to converge within MAXIT iterations.
!C
      ITER = MAXIT
      IF (.NOT. LOOP2) THEN
        LOOP2 = .TRUE.
        GO TO 2
      ENDIF
!C
!C Store convergence parameters.
!C
    7 WK(1) = DYP
      WK(2) = DSMAX
      IF (LOOP2) ITER = ITER + MAXIT
!C
!C No error encountered.
!C
   10 IER = ITER
      RETURN
!C
!C Invalid input parameter N, NCD, or IENDC.
!C
   11 IER = -1
      RETURN
!C
!C LWK too small.
!C
   12 IER = -2
      RETURN
!C
!C Abscissae are not strictly increasing.
!C
   14 IER = -4
      RETURN
      END SUBROUTINE
      
      SUBROUTINE TSVAL1 (N,X,Y,YP,SIGMA,IFLAG,NE,TE, V,                 &
     &                   IER)
      INTEGER N, IFLAG, NE, IER
      DOUBLE PRECISION X(N), Y(N), YP(N), SIGMA(N), TE(NE),             &
     &                 V(NE)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This subroutine evaluates a Hermite interpolatory ten-
!C sion spline H or its first or second derivative at a set
!C of points TE.
!C
!C   Note that a large tension factor in SIGMA may cause
!C underflow.  The result is assumed to be zero.  If not the
!C default, this may be specified by either a compiler option
!C or operating system option.
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing the abscissae.
!C           These must be in strictly increasing order:
!C           X(I) < X(I+1) for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values or
!C           function values returned by Subroutine SMCRV.
!C           Y(I) = H(X(I)) for I = 1,...,N.
!C
!C       YP = Array of length N containing first deriva-
!C            tives.  YP(I) = HP(X(I)) for I = 1,...,N, where
!C            HP denotes the derivative of H.
!C
!C       SIGMA = Array of length N-1 containing tension fac-
!C               tors whose absolute values determine the
!C               balance between cubic and linear in each
!C               interval.  SIGMA(I) is associated with int-
!C               erval (I,I+1) for I = 1,...,N-1.
!C
!C       IFLAG = Output option indicator:
!C               IFLAG = 0 if values of H are to be computed.
!C               IFLAG = 1 if first derivative values are to
!C                         be computed.
!C               IFLAG = 2 if second derivative values are to
!C                         be computed.
!C
!C       NE = Number of evaluation points.  NE > 0.
!C
!C       TE = Array of length NE containing the evaluation
!C            points.  The sequence should be strictly in-
!C            creasing for maximum efficiency.  Extrapolation
!C            is performed if a point is not in the interval
!C            [X(1),X(N)].
!C
!C The above parameters are not altered by this routine.
!C
!C       V = Array of length at least NE.
!C
!C On output:
!C
!C       V = Array of function, first derivative, or second
!C           derivative values at the evaluation points un-
!C           less IER < 0.  If IER = -1, V is not altered.
!C           If IER = -2, V may be only partially defined.
!C
!C       IER = Error indicator:
!C             IER = 0  if no errors were encountered and
!C                      no extrapolation occurred.
!C             IER > 0  if no errors were encountered but
!C                      extrapolation was required at IER
!C                      points.
!C             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or
!C                      NE < 1.
!C             IER = -2 if the abscissae are not in strictly
!C                      increasing order.  (This error will
!C                      not necessarily be detected.)
!C
!C Modules required by TSVAL1:  HPPVAL, HPVAL, HVAL, INTRVL,
!C                                SNHCSH
!C
!C***********************************************************
!C
      INTEGER I, IERR, IFLG, NVAL, NX
! in a module these need not to be declared here      DOUBLE PRECISION HPPVAL, HPVAL, HVAL
!C
      IFLG = IFLAG
      NVAL = NE
!C
!C Test for invalid input.
!C
      IF (N .LT. 2  .OR.  IFLG .LT. 0  .OR.  IFLG .GT. 2                &
     &    .OR.  NVAL .LT. 1) GO TO 2
!C
!C Initialize the number of extrapolation points NX and
!C   loop on evaluation points.
!C
      NX = 0
      DO 1 I = 1,NVAL
        IF (IFLG .EQ. 0) THEN
          V(I) = HVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ELSEIF (IFLG .EQ. 1) THEN
          V(I) = HPVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ELSE
          V(I) = HPPVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ENDIF
        IF (IERR .LT. 0) GO TO 3
        NX = NX + IERR
    1   CONTINUE
!C
!C No errors encountered.
!C
      IER = NX
      RETURN
!C
!C N, IFLAG, or NE is outside its valid range.
!C
    2 IER = -1
      RETURN
!C
!C X is not strictly increasing.
!C
    3 IER = -2
      RETURN
      END SUBROUTINE
      
      SUBROUTINE YPC1 (N,X,Y, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This subroutine employs a three-point quadratic interpo-
!C lation method to compute local derivative estimates YP
!C associated with a set of data points.  The interpolation
!C formula is the monotonicity-constrained parabolic method
!C described in the reference cited below.  A Hermite int-
!C erpolant of the data values and derivative estimates pre-
!C serves monotonicity of the data.  Linear interpolation is
!C used if N = 2.  The method is invariant under a linear
!C scaling of the coordinates but is not additive.
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing a strictly in-
!C           creasing sequence of abscissae:  X(I) < X(I+1)
!C           for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values asso-
!C           ciated with the abscissae.
!C
!C Input parameters are not altered by this routine.
!C
!C On output:
!C
!C       YP = Array of length N containing estimated deriv-
!C            atives at the abscissae unless IER .NE. 0.
!C            YP is not altered if IER = 1, and is only par-
!C            tially defined if IER > 1.
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered.
!C             IER = 1 if N < 2.
!C             IER = I if X(I) .LE. X(I-1) for some I in the
!C                     range 2,...,N.
!C
!C Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
!C               Cubic Interpolation",  LA-8796-MS, Los
!C               Alamos National Lab, Feb. 1982.
!C
!C Modules required by YPC1:  None
!C
!C Intrinsic functions called by YPC1:  ABS, MAX, MIN, SIGN
!C
!C***********************************************************
!C
      INTEGER I, NM1
      DOUBLE PRECISION ASI, ASIM1, DX2, DXI, DXIM1, S2, SGN,            &
     &                 SI, SIM1, T
!C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 2
      I = 1
      DXI = X(2) - X(1)
      IF (DXI .LE. 0.D0) GO TO 3
      SI = (Y(2)-Y(1))/DXI
      IF (NM1 .EQ. 1) THEN
!C
!C Use linear interpolation for N = 2.
!C
        YP(1) = SI
        YP(2) = SI
        IER = 0
        RETURN
      ENDIF
!C
!C N .GE. 3.  YP(1) = S1 + DX1*(S1-S2)/(DX1+DX2) unless this
!C   results in YP(1)*S1 .LE. 0 or abs(YP(1)) > 3*abs(S1).
!C
      I = 2
      DX2 = X(3) - X(2)
      IF (DX2 .LE. 0.D0) GO TO 3
      S2 = (Y(3)-Y(2))/DX2
      T = SI + DXI*(SI-S2)/(DXI+DX2)
      IF (SI .GE. 0.D0) THEN
        YP(1) = MIN(MAX(0.D0,T), 3.D0*SI)
      ELSE
        YP(1) = MAX(MIN(0.D0,T), 3.D0*SI)
      ENDIF
!C
!C YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the
!C   constraint that YP(I) has the sign of either SIM1 or
!C   SI, whichever has larger magnitude, and abs(YP(I)) .LE.
!C   3*min(abs(SIM1),abs(SI)).
!C
      DO 1 I = 2,NM1
        DXIM1 = DXI
        DXI = X(I+1) - X(I)
        IF (DXI .LE. 0.D0) GO TO 3
        SIM1 = SI
        SI = (Y(I+1)-Y(I))/DXI
        T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI)
        ASIM1 = DABS(SIM1)
        ASI = DABS(SI)
        SGN = SIGN(1.D0,SI)
        IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
        IF (SGN .GT. 0.D0) THEN
          YP(I) = MIN(MAX(0.D0,T),3.D0*MIN(ASIM1,ASI))
        ELSE
          YP(I) = MAX(MIN(0.D0,T),-3.D0*MIN(ASIM1,ASI))
        ENDIF
    1   CONTINUE
!C
!C YP(N) = SNM1 + DXNM1*(SNM1-SNM2)/(DXNM2+DXNM1) subject to
!C   the constraint that YP(N) has the sign of SNM1 and
!C   abs(YP(N)) .LE. 3*abs(SNM1).  Note that DXI = DXNM1 and
!C   SI = SNM1.
!C
      T = SI + DXI*(SI-SIM1)/(DXIM1+DXI)
      IF (SI .GE. 0.D0) THEN
        YP(N) = MIN(MAX(0.D0,T), 3.D0*SI)
      ELSE
        YP(N) = MAX(MIN(0.D0,T), 3.D0*SI)
      ENDIF
      IER = 0
      RETURN
!C
!C N is outside its valid range.
!C
    2 IER = 1
      RETURN
!C
!C X(I+1) .LE. X(I).
!C
    3 IER = I + 1
      RETURN
      END SUBROUTINE
      
      SUBROUTINE YPC1P (N,X,Y, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This subroutine employs a three-point quadratic interpo-
!C lation method to compute local derivative estimates YP
!C associated with a set of N data points (X(I),Y(I)).  It
!C is assumed that Y(N) = Y(1), and YP(N) = YP(1) on output.
!C Thus, a Hermite interpolant H(x) defined by the data
!C points and derivative estimates is periodic with period
!C X(N)-X(1).  The derivative-estimation formula is the
!C monotonicity-constrained parabolic fit described in the
!C reference cited below:  H(x) is monotonic in intervals in
!C which the data is monotonic.  The method is invariant
!C under a linear scaling of the coordinates but is not
!C additive.
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 3.
!C
!C       X = Array of length N containing a strictly in-
!C           creasing sequence of abscissae:  X(I) < X(I+1)
!C           for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values asso-
!C           ciated with the abscissae.  Y(N) is set to Y(1)
!C           on output unless IER = 1.
!C
!C   Input parameters, other than Y(N), are not altered by
!C this routine.
!C
!C On output:
!C
!C       YP = Array of length N containing estimated deriv-
!C            atives at the abscissae unless IER .NE. 0.
!C            YP is not altered if IER = 1, and is only par-
!C            tially defined if IER > 1.
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered.
!C             IER = 1 if N < 3.
!C             IER = I if X(I) .LE. X(I-1) for some I in the
!C                     range 2,...,N.
!C
!C Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
!C               Cubic Interpolation",  LA-8796-MS, Los
!C               Alamos National Lab, Feb. 1982.
!C
!C Modules required by YPC1P:  None
!C
!C Intrinsic functions called by YPC1P:  ABS, MAX, MIN, SIGN
!C
!C***********************************************************
!C
      INTEGER I, NM1
      DOUBLE PRECISION ASI, ASIM1, DXI, DXIM1, SGN, SI,                 &
     &                 SIM1, T
!C
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 2
      Y(N) = Y(1)
!C
!C Initialize for loop on interior points.
!C
      I = 1
      DXI = X(2) - X(1)
      IF (DXI .LE. 0.D0) GO TO 3
      SI = (Y(2)-Y(1))/DXI
!C
!C YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the
!C   constraint that YP(I) has the sign of either SIM1 or
!C   SI, whichever has larger magnitude, and abs(YP(I)) .LE.
!C   3*min(abs(SIM1),abs(SI)).
!C
      DO 1 I = 2,NM1
        DXIM1 = DXI
        DXI = X(I+1) - X(I)
        IF (DXI .LE. 0.D0) GO TO 3
        SIM1 = SI
        SI = (Y(I+1)-Y(I))/DXI
        T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI)
        ASIM1 = DABS(SIM1)
        ASI = DABS(SI)
        SGN = SIGN(1.D0,SI)
        IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
        IF (SGN .GT. 0.D0) THEN
          YP(I) = MIN(MAX(0.D0,T),3.D0*MIN(ASIM1,ASI))
        ELSE
          YP(I) = MAX(MIN(0.D0,T),-3.D0*MIN(ASIM1,ASI))
        ENDIF
    1   CONTINUE
!C
!C YP(N) = YP(1), I = 1, and IM1 = N-1.
!C
      DXIM1 = DXI
      DXI = X(2) - X(1)
      SIM1 = SI
      SI = (Y(2) - Y(1))/DXI
      T = (DXIM1*SI + DXI*SIM1)/(DXIM1+DXI)
      ASIM1 = DABS(SIM1)
      ASI = DABS(SI)
      SGN = SIGN(1.D0,SI)
      IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
      IF (SGN .GT. 0.D0) THEN
        YP(1) = MIN(MAX(0.D0,T), 3.D0*MIN(ASIM1,ASI))
      ELSE
        YP(1) = MAX(MIN(0.D0,T), -3.D0*MIN(ASIM1,ASI))
      ENDIF
      YP(N) = YP(1)
!C
!C No error encountered.
!C
      IER = 0
      RETURN
!C
!C N is outside its valid range.
!C
    2 IER = 1
      RETURN
!C
!C X(I+1) .LE. X(I).
!C
    3 IER = I + 1
      RETURN
      END SUBROUTINE
      
      SUBROUTINE YPC2 (N,X,Y,SIGMA,ISL1,ISLN,BV1,BVN,                   &
     &                 WK, YP,IER)
      INTEGER N, ISL1, ISLN, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), BV1, BVN,                  &
     &                 WK(N), YP(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This subroutine solves a linear system for a set of
!C first derivatives YP associated with a Hermite interpola-
!C tory tension spline H(x).  The derivatives are chosen so
!C that H(x) has two continuous derivatives for all x and H
!C satisfies user-specified end conditions.
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing a strictly in-
!C           creasing sequence of abscissae:  X(I) < X(I+1)
!C           for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values asso-
!C           ciated with the abscissae.  H(X(I)) = Y(I) for
!C           I = 1,...,N.
!C
!C       SIGMA = Array of length N-1 containing tension
!C               factors.  SIGMA(I) is associated with inter-
!C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
!C               SIGMA(I) = 0, H is the Hermite cubic interp-
!C               olant of the data values and computed deriv-
!C               atives at X(I) and X(I+1), and if all
!C               tension factors are zero, H is the !C-2 cubic
!C               spline interpolant which satisfies the end
!C               conditions.
!C
!C       ISL1 = Option indicator for the condition at X(1):
!C              ISL1 = 0 if YP(1) is to be estimated inter-
!C                       nally by a constrained parabolic
!C                       fit to the first three points.
!C                       This is identical to the method used
!C                       by Subroutine YPC1.  BV1 is not used
!C                       in this case.
!C              ISL1 = 1 if the first derivative of H at X(1)
!C                       is specified by BV1.
!C              ISL1 = 2 if the second derivative of H at
!C                       X(1) is specified by BV1.
!C              ISL1 = 3 if YP(1) is to be estimated inter-
!C                       nally from the derivative of the
!C                       tension spline (using SIGMA(1))
!C                       which interpolates the first three
!C                       data points and has third derivative
!C                       equal to zero at X(1).  Refer to
!C                       ENDSLP.  BV1 is not used in this
!C                       case.
!C
!C       ISLN = Option indicator for the condition at X(N):
!C              ISLN = 0 if YP(N) is to be estimated inter-
!C                       nally by a constrained parabolic
!C                       fit to the last three data points.
!C                       This is identical to the method used
!C                       by Subroutine YPC1.  BVN is not used
!C                       in this case.
!C              ISLN = 1 if the first derivative of H at X(N)
!C                       is specified by BVN.
!C              ISLN = 2 if the second derivative of H at
!C                       X(N) is specified by BVN.
!C              ISLN = 3 if YP(N) is to be estimated inter-
!C                       nally from the derivative of the
!C                       tension spline (using SIGMA(N-1))
!C                       which interpolates the last three
!C                       data points and has third derivative
!C                       equal to zero at X(N).  Refer to
!C                       ENDSLP.  BVN is not used in this
!C                       case.
!C
!C       BV1,BVN = Boundary values or dummy parameters as
!C                 defined by ISL1 and ISLN.
!C
!C The above parameters are not altered by this routine.
!C
!C       WK = Array of length at least N-1 to be used as
!C            temporary work space.
!C
!C       YP = Array of length .GE. N.
!C
!C On output:
!C
!C       YP = Array containing derivatives of H at the
!C            abscissae.  YP is not defined if IER .NE. 0.
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered.
!C             IER = 1 if N, ISL1, or ISLN is outside its
!C                     valid range.
!C             IER = I if X(I) .LE. X(I-1) for some I in the
!C                     range 2,...,N.
!C
!C Modules required by YPC2:  ENDSLP, SNHCSH, YPCOEF
!C
!C Intrinsic function called by YPC2:  ABS
!C
!C***********************************************************
!C
      INTEGER I, NM1, NN
      DOUBLE PRECISION D, D1, D2, DX, R1, R2, S, SD1, SD2,              &
     &                 SIG, YP1, YPN
! in a module these need not to be declared here      DOUBLE PRECISION ENDSLP
!C
      NN = N
      IF (NN .LT. 2  .OR.  ISL1 .LT. 0  .OR.  ISL1 .GT. 3               &
     &    .OR.  ISLN .LT. 0  .OR.  ISLN .GT. 3) GO TO 3
      NM1 = NN - 1
!C
!C Set YP1 and YPN to the endpoint values.
!C
      IF (ISL1 .EQ. 0) THEN
        IF (NN .GT. 2) YP1 = ENDSLP (X(1),X(2),X(3),Y(1),               &
     &                               Y(2),Y(3),0.D0)
      ELSEIF (ISL1 .NE. 3) THEN
        YP1 = BV1
      ELSE
        IF (NN .GT. 2) YP1 = ENDSLP (X(1),X(2),X(3),Y(1),               &
     &                               Y(2),Y(3),SIGMA(1))
      ENDIF
      IF (ISLN .EQ. 0) THEN
        IF (NN .GT. 2) YPN = ENDSLP (X(NN),X(NM1),X(NN-2),              &
     &                            Y(NN),Y(NM1),Y(NN-2),0.D0)
      ELSEIF (ISLN .NE. 3) THEN
        YPN = BVN
      ELSE
        IF (NN .GT. 2) YPN = ENDSLP (X(NN),X(NM1),X(NN-2),              &
     &                      Y(NN),Y(NM1),Y(NN-2),SIGMA(NM1))
      ENDIF
!C
!C Solve the symmetric positive-definite tridiagonal linear
!C   system.  The forward elimination step consists of div-
!C   iding each row by its diagonal entry, then introducing a
!C   zero below the diagonal.  This requires saving only the
!C   superdiagonal (in WK) and the right hand side (in YP).
!C
      I = 1
      DX = X(2) - X(1)
      IF (DX .LE. 0.D0) GO TO 4
      S = (Y(2)-Y(1))/DX
      IF (NN .EQ. 2) THEN
        IF (ISL1 .EQ. 0  .OR.  ISL1 .EQ. 3) YP1 = S
        IF (ISLN .EQ. 0  .OR.  ISLN .EQ. 3) YPN = S
      ENDIF
!C
!C Begin forward elimination.
!C
      SIG = DABS(SIGMA(1))
      CALL YPCOEF (SIG,DX, D1,SD1)
      R1 = (SD1+D1)*S
      WK(1) = 0.D0
      YP(1) = YP1
      IF (ISL1 .EQ. 2) THEN
        WK(1) = SD1/D1
        YP(1) = (R1-YP1)/D1
      ENDIF
      DO 1 I = 2,NM1
        DX = X(I+1) - X(I)
        IF (DX .LE. 0.D0) GO TO 4
        S = (Y(I+1)-Y(I))/DX
        SIG = DABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D2,SD2)
        R2 = (SD2+D2)*S
        D = D1 + D2 - SD1*WK(I-1)
        WK(I) = SD2/D
        YP(I) = (R1 + R2 - SD1*YP(I-1))/D
        D1 = D2
        SD1 = SD2
        R1 = R2
    1   CONTINUE
      D = D1 - SD1*WK(NM1)
      YP(NN) = YPN
      IF (ISLN .EQ. 2) YP(NN) = (R1 + YPN - SD1*YP(NM1))/D
!C
!C Back substitution:
!C
      DO 2 I = NM1,1,-1
        YP(I) = YP(I) - WK(I)*YP(I+1)
    2   CONTINUE
      IER = 0
      RETURN
!C
!C Invalid integer input parameter.
!C
    3 IER = 1
      RETURN
!C
!C Abscissae out of order or duplicate points.
!C
    4 IER = I + 1
      RETURN
      END SUBROUTINE
      
      SUBROUTINE YPC2P (N,X,Y,SIGMA,WK, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), WK(*), YP(N)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   06/10/92
!C
!C   This subroutine solves a linear system for a set of
!C first derivatives YP associated with a Hermite interpola-
!C tory tension spline H(x).  The derivatives are chosen so
!C that H(x) has two continuous derivatives for all x, and H
!C satisfies periodic end conditions:  first and second der-
!C ivatives of H at X(1) agree with those at X(N), and thus
!C the length of a period is X(N) - X(1).  It is assumed that
!C Y(N) = Y(1), and Y(N) is not referenced.
!C
!C On input:
!C
!C       N = Number of data points.  N .GE. 3.
!C
!C       X = Array of length N containing a strictly in-
!C           creasing sequence of abscissae:  X(I) < X(I+1)
!C           for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values asso-
!C           ciated with the abscissae.  H(X(I)) = Y(I) for
!C           I = 1,...,N.
!C
!C       SIGMA = Array of length N-1 containing tension
!C               factors.  SIGMA(I) is associated with inter-
!C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
!C               SIGMA(I) = 0, H is the Hermite cubic interp-
!C               olant of the data values and computed deriv-
!C               atives at X(I) and X(I+1), and if all
!C               tension factors are zero, H is the !C-2 cubic
!C               spline interpolant which satisfies the end
!C               conditions.
!C
!C The above parameters are not altered by this routine.
!C
!C       WK = Array of length at least 2N-2 to be used as
!C            temporary work space.
!C
!C       YP = Array of length .GE. N.
!C
!C On output:
!C
!C       YP = Array containing derivatives of H at the
!C            abscissae.  YP is not defined if IER .NE. 0.
!C
!C       IER = Error indicator:
!C             IER = 0 if no errors were encountered.
!C             IER = 1 if N is outside its valid range.
!C             IER = I if X(I) .LE. X(I-1) for some I in the
!C                     range 2,...,N.
!C
!C Modules required by YPC2P:  SNHCSH, YPCOEF
!C
!C Intrinsic function called by YPC2P:  ABS
!C
!C***********************************************************
!C
      INTEGER I, NM1, NM2, NM3, NN, NP1, NPI
      DOUBLE PRECISION D, D1, D2, DIN, DNM1, DX, R1, R2,                &
     &                 RNM1, S, SD1, SD2, SDNM1, SIG, YPNM1
!C
      NN = N
      IF (NN .LT. 3) GO TO 4
      NM1 = NN - 1
      NM2 = NN - 2
      NM3 = NN - 3
      NP1 = NN + 1
!C
!C The system is order N-1, symmetric, positive-definite, and
!C   tridiagonal except for nonzero elements in the upper
!C   right and lower left corners.  The forward elimination
!C   step zeros the subdiagonal and divides each row by its
!C   diagonal entry for the first N-2 rows.  The superdiago-
!C   nal is stored in WK(I), the negative of the last column
!C   (fill-in) in WK(N+I), and the right hand side in YP(I)
!C   for I = 1,...,N-2.
!C
      I = NM1
      DX = X(NN) - X(NM1)
      IF (DX .LE. 0.D0) GO TO 5
      S = (Y(1)-Y(NM1))/DX
      SIG = DABS(SIGMA(NM1))
      CALL YPCOEF (SIG,DX, DNM1,SDNM1)
      RNM1 = (SDNM1+DNM1)*S
      I = 1
      DX = X(2) - X(1)
      IF (DX .LE. 0.D0) GO TO 5
      S = (Y(2)-Y(1))/DX
      SIG = DABS(SIGMA(1))
      CALL YPCOEF (SIG,DX, D1,SD1)
      R1 = (SD1+D1)*S
      D = DNM1 + D1
      WK(1) = SD1/D
      WK(NP1) = -SDNM1/D
      YP(1) = (RNM1+R1)/D
      DO 1 I = 2,NM2
        DX = X(I+1) - X(I)
        IF (DX .LE. 0.D0) GO TO 5
        S = (Y(I+1)-Y(I))/DX
        SIG = DABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D2,SD2)
        R2 = (SD2+D2)*S
        D = D1 + D2 - SD1*WK(I-1)
        DIN = 1.D0/D
        WK(I) = SD2*DIN
        NPI = NN + I
        WK(NPI) = -SD1*WK(NPI-1)*DIN
        YP(I) = (R1 + R2 - SD1*YP(I-1))*DIN
        SD1 = SD2
        D1 = D2
        R1 = R2
    1   CONTINUE
!C
!C The backward elimination step zeros the superdiagonal
!C   (first N-3 elements).  WK(I) and YP(I) are overwritten
!C   with the negative of the last column and the new right
!C   hand side, respectively, for I = N-2, N-3, ..., 1.
!C
      NPI = NN + NM2
      WK(NM2) = WK(NPI) - WK(NM2)
      DO 2 I = NM3,1,-1
        YP(I) = YP(I) - WK(I)*YP(I+1)
        NPI = NN + I
        WK(I) = WK(NPI) - WK(I)*WK(I+1)
    2   CONTINUE
!C
!C Solve the last equation for YP(N-1).
!C
      YPNM1 = (R1 + RNM1 - SDNM1*YP(1) - SD1*YP(NM2))/                  &
     &        (D1 + DNM1 + SDNM1*WK(1) + SD1*WK(NM2))
!C
!C Back substitute for the remainder of the solution
!C   components.
!C
      YP(NM1) = YPNM1
      DO 3 I = 1,NM2
        YP(I) = YP(I) + WK(I)*YPNM1
    3   CONTINUE
!C
!C YP(N) = YP(1).
!C
      YP(N) = YP(1)
      IER = 0
      RETURN
!C
!C N is outside its valid range.
!C
    4 IER = 1
      RETURN
!C
!C Abscissae out of order or duplicate points.
!C
    5 IER = I + 1
      RETURN
      END SUBROUTINE

      SUBROUTINE YPCOEF (SIGMA,DX, D,SD)
      DOUBLE PRECISION SIGMA, DX, D, SD
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/17/96
!C
!C   This subroutine computes the coefficients of the deriva-
!C tives in the symmetric diagonally dominant tridiagonal
!C system associated with the !C-2 derivative estimation pro-
!C cedure for a Hermite interpolatory tension spline.
!C
!C On input:
!C
!C       SIGMA = Nonnegative tension factor associated with
!C               an interval.
!C
!C       DX = Positive interval width.
!C
!C Input parameters are not altered by this routine.
!C
!C On output:
!C
!C       D = Component of the diagonal term associated with
!C           the interval.  D = SIG*(SIG*COSHM(SIG) -
!C           SINHM(SIG))/(DX*E), where SIG = SIGMA and E =
!C           SIG*SINH(SIG) - 2*COSHM(SIG).
!C
!C       SD = Subdiagonal (superdiagonal) term.  SD = SIG*
!C            SINHM(SIG)/E.
!C
!C Module required by YPCOEF:  SNHCSH
!C
!C Intrinsic function called by YPCOEF:  EXP
!C
!C***********************************************************
!C
      DOUBLE PRECISION COSHM, COSHMM, E, EMS, SCM, SIG,                 &
     &                 SINHM, SSINH, SSM
!C
      SIG = SIGMA
      IF (SIG .LT. 1.D-9) THEN
!C
!C SIG = 0:  cubic interpolant.
!C
        D = 4.D0/DX
        SD = 2.D0/DX
      ELSEIF (SIG .LE. .5D0) THEN
!C
!C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!C                      cancellation error in the hyperbolic
!C                      functions when SIGMA is small.
!C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = (SIG*SINHM - COSHMM - COSHMM)*DX
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
!C
!C SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
!C            to avoid overflow when SIGMA is large.
!C
        EMS = EXP(-SIG)
        SSINH = 1.D0 - EMS*EMS
        SSM = SSINH - 2.D0*SIG*EMS
        SCM = (1.D0-EMS)*(1.D0-EMS)
        E = (SIG*SSINH - SCM - SCM)*DX
        D = SIG*(SIG*SCM-SSM)/E
        SD = SIG*SSM/E
      ENDIF
      RETURN
      END SUBROUTINE

      DOUBLE PRECISION FUNCTION my_HPVAL (T,X,Y,YP,SIGMA,IER)
      INTEGER IER
      DOUBLE PRECISION T, X(2), Y(2), YP(2), SIGMA(2)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu
!C                                                   11/17/96
!C
!C   This function evaluates the first derivative HP of a
!C Hermite interpolatory tension spline H at a point T.
!C
!C On input:
!C
!C       T = Point at which HP is to be evaluated.  Extrapo-
!C           lation is performed if T < X(1) or T > X(N).
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing the abscissae.
!C           These must be in strictly increasing order:
!C           X(I) < X(I+1) for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values.
!C           H(X(I)) = Y(I) for I = 1,...,N.
!C
!C       YP = Array of length N containing first deriva-
!C            tives.  HP(X(I)) = YP(I) for I = 1,...,N.
!C
!C       SIGMA = Array of length N-1 containing tension fac-
!C               tors whose absolute values determine the
!C               balance between cubic and linear in each
!C               interval.  SIGMA(I) is associated with int-
!C               erval (I,I+1) for I = 1,...,N-1.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0  if no errors were encountered and
!C                      X(1) .LE. T .LE. X(N).
!C             IER = 1  if no errors were encountered and
!C                      extrapolation was necessary.
!C             IER = -1 if N < 2.
!C             IER = -2 if the abscissae are not in strictly
!C                      increasing order.  (This error will
!C                      not necessarily be detected.)
!C
!C       HPVAL = Derivative value HP(T), or zero if IER < 0.
!C
!C Modules required by HPVAL:  INTRVL, SNHCSH
!C
!C Intrinsic functions called by HPVAL:  ABS, EXP
!C
!C***********************************************************
!C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,             &
     &                 DX, E, E1, E2, EMS, S, S1, SB1, SB2,             &
     &                 SBIG, SIG, SINH2, SM, SM2, TM
! in a module these need not to be declared here      INTEGER INTRVL
!C
      DATA SBIG/85.D0/
!      IF (N .LT. 2) GO TO 1
      IER = 0
      I = 1
      IP1 = I + 1
!C
!C Compute interval width DX, local coordinates B1 and B2,
!C   and second differences D1 and D2.
!C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      B1 = (X(IP1) - T)/DX
      B2 = 1.D0 - B1
      S1 = YP(I)
      S = (Y(IP1)-Y(I))/DX
      D1 = S - S1
      D2 = YP(IP1) - S
      SIG = DABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
!C
!C SIG = 0:  H is the Hermite cubic interpolant.
!C
        my_HPVAL = S1 + B2*(D1 + D2 - 3.D0*B1*(D2-D1))
      ELSEIF (SIG .LE. .5D0) THEN
!C
!C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!C   cancellation error in the hyperbolic functions.
!C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH2 = SM2 + SB2
        E = SIG*SM - CMM - CMM
        my_HPVAL = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) +                    &
     &                SIG*(CM*SINH2-(SM+SIG)*CM2)*D1)/E
      ELSE
!C
!C SIG > .5:  use negative exponentials in order to avoid
!C   overflow.  Note that EMS = EXP(-SIG).  In the case of
!C   extrapolation (negative B1 or B2), H is approximated by
!C   a linear function if -SIG*B1 or -SIG*B2 is large.
!C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          my_HPVAL = S
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          E = TM*(SIG*(1.D0+EMS) - TM - TM)
          my_HPVAL = S + (TM*((E2-E1)*(D1+D2) + TM*(D1-D2)) +           &
     &            SIG*((E1*EMS-E2)*D1 + (E1-E2*EMS)*D2))/E
        ENDIF
      ENDIF
      RETURN
!C
!C N is outside its valid range.
!C
!    1 my_HPVAL = 0.D0
!      IER = -1
!      RETURN
!C
!C X(I) .GE. X(I+1).
!C
    2 my_HPVAL = 0.D0
      IER = -2
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION my_HVAL (T,X,Y,YP,SIGMA,IER)
      INTEGER IER
      DOUBLE PRECISION T, X(2), Y(2), YP(2), SIGMA(2)
!C
!C***********************************************************
!C
!C                                                From TSPACK
!C                                            Robert J. Renka
!C                                  Dept. of Computer Science
!C                                       Univ. of North Texas
!C                                           renka@cs.unt.edu

!C                                                   11/17/96
!C
!C   This function evaluates a Hermite interpolatory tension
!C spline H at a point T.  Note that a large value of SIGMA
!C may cause underflow.  The result is assumed to be zero.
!C
!C   Given arrays X, Y, YP, and SIGMA of length NN, if T is
!C known to lie in the interval (X(I),X(J)) for some I < J,
!C a gain in efficiency can be achieved by calling this
!C function with N = J+1-I (rather than NN) and the I-th
!C components of the arrays (rather than the first) as par-
!C ameters.
!C
!C On input:
!C
!C       T = Point at which H is to be evaluated.  Extrapo-
!C           lation is performed if T < X(1) or T > X(N).
!C
!C       N = Number of data points.  N .GE. 2.
!C
!C       X = Array of length N containing the abscissae.
!C           These must be in strictly increasing order:
!C           X(I) < X(I+1) for I = 1,...,N-1.
!C
!C       Y = Array of length N containing data values.
!C           H(X(I)) = Y(I) for I = 1,...,N.
!C
!C       YP = Array of length N containing first deriva-
!C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
!C            HP denotes the derivative of H.
!C
!C       SIGMA = Array of length N-1 containing tension fac-
!C               tors whose absolute values determine the
!C               balance between cubic and linear in each
!C               interval.  SIGMA(I) is associated with int-
!C               erval (I,I+1) for I = 1,...,N-1.
!C
!C Input parameters are not altered by this function.
!C
!C On output:
!C
!C       IER = Error indicator:
!C             IER = 0  if no errors were encountered and
!C                      X(1) .LE. T .LE. X(N).
!C             IER = 1  if no errors were encountered and
!C                      extrapolation was necessary.
!C             IER = -1 if N < 2.
!C             IER = -2 if the abscissae are not in strictly
!C                      increasing order.  (This error will
!C                      not necessarily be detected.)
!C
!C       HVAL = Function value H(T), or zero if IER < 0.
!C
!C Modules required by HVAL:  INTRVL, SNHCSH
!C
!C Intrinsic functions called by HVAL:  ABS, EXP
!C
!C***********************************************************
!C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,             &
     &                 DX, E, E1, E2, EMS, S, S1, SB1, SB2,             &
     &                 SBIG, SIG, SM, SM2, TM, TP, TS, U, Y1
! in a module these need not to be declared here      INTEGER INTRVL
!C
      DATA SBIG/85.D0/
!      IF (N .LT. 2) GO TO 1
      IER = 0
      I = 1
      IP1 = I + 1
!C
!C Compute interval width DX, local coordinates B1 and B2,
!C   and second differences D1 and D2.
!C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      U = T - X(I)
      B2 = U/DX
      B1 = 1.D0 - B2
      Y1 = Y(I)
      S1 = YP(I)
      S = (Y(IP1)-Y1)/DX
      D1 = S - S1
      D2 = YP(IP1) - S
      SIG = DABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
!C
!C SIG = 0:  H is the Hermite cubic interpolant.
!C
        my_HVAL = Y1 + U*(S1 + B2*(D1 + B1*(D1-D2)))
      ELSEIF (SIG .LE. .5D0) THEN
!C
!C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!C   cancellation error in the hyperbolic functions.
!C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        my_HVAL = Y1 + S1*U + DX*((CM*SM2-SM*CM2)*(D1+D2) +             &
     &                         SIG*(CM*CM2-(SM+SIG)*SM2)*D1)            &
     &                         /(SIG*E)
      ELSE
!C
!C SIG > .5:  use negative exponentials in order to avoid
!C   overflow.  Note that EMS = EXP(-SIG).  In the case of
!C   extrapolation (negative B1 or B2), H is approximated by
!C   a linear function if -SIG*B1 or -SIG*B2 is large.
!C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          my_HVAL = Y1 + S*U
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TS = TM*TM
          TP = 1.D0 + EMS
          E = TM*(SIG*TP - TM - TM)
          my_HVAL = Y1 + S*U + DX*(TM*(TP-E1-E2)*(D1+D2) + SIG*         &
     &                         ((E2+EMS*(E1-2.D0)-B1*TS)*D1+            &
     &                        (E1+EMS*(E2-2.D0)-B2*TS)*D2))/            &
     &                        (SIG*E)
        ENDIF
      ENDIF
      RETURN
!C
!C N is outside its valid range.
!C
!    1 my_HVAL = 0.D0
!      IER = -1
!      RETURN
!C
!C X(I) .GE. X(I+1).
!C
    2 my_HVAL = 0.D0
      IER = -2
      RETURN
      END FUNCTION

end module TSPACK

      subroutine EvalTabulatedFunction(inverse,n,ind1,ind2,ind3,node,   &
     &                                  sptab,ientrytab,xe,ye,dyedxe, iWhat)
      use doln
      use doTSPACK
      use TSPACK
      implicit none
      include 'arrays.fi'

      integer k, klo, khi, n, inverse, maxtry, ntry, ind1, ind2, ind3, iWhat
      real(8) xe, ye, dyedxe
      real(8) x1, x2, f1, f2, d1, d2, h, delta, del1, del2, c2, c2t2,   &
     &        c3, c3t3, xx
      real(8) xlow, xhigh, xtry, ytry, xe_local
      real(8) :: sptab(7,macp,matab)
      real(8) :: x(2), y(2), yp(2), sig(2)
      integer ientrytab(macp,0:matabentries), node    
      save

!     'ye' is only used as 'inverse' = 1. The option that requires this is not yet operationel
!     (combination of sophys-tables (swsophy=1) and swdiscrvert=1).
! Therefore, 'ye' initial gets a dummy value for Forecheck
!      if (inverse == 1) ye = -99.d0
!      if (inverse == 0) xe = -99.d0
       ye = 0.d0

!      sptab(1,node,i):    h
!      sptab(2,node,i):    Theta
!      sptab(3,node,i):    K
!      sptab(4,node,i):    dTheta / dh
!      sptab(5,node,i):    dK / dh
!      sptab(6,node,i):    sigma (Theta(h))
!      sptab(7,node,i):    sigma (K(h))

!     bi-section search method for finding the array entry

      if(inverse .eq. 0)then
         if(xe.ge.-1.0d-9) then
            stop 'Invalid call of Function EvalTabulatedFunction'
         end if
         k   = int(1000*(log10(-xe)+1.d0))+4001
         if(k.lt.1) k=1
         klo = ientrytab(node,k) 
         khi = klo + 1

         if (do_ln_trans) then
            xe_local = -(dlog(-xe+1.0d0))
         else
            xe_local = xe
         end if
         if (xe < -16.5d0) then
            continue
         end if
!     bounds of interval
         x1 = sptab(ind1,node,klo)
         x2 = sptab(ind1,node,khi)

!         write (234,'(I5,4F20.10)') klo, xe, xe_local, x1, x2

         
         if( (xe_local-x1) .lt. 0.0d0 )then
            klo = klo -1
            khi = khi -1
            x1 = sptab(ind1,node,klo)
            x2 = sptab(ind1,node,khi)
         else if( (xe_local-x1) .ge. (x2-x1) )then
            klo = klo +1
            khi = khi +1
            x1 = sptab(ind1,node,klo)
            x2 = sptab(ind1,node,khi)
         end if

         if (use_TSPACK) then
            x(1:2)   = sptab(ind1,node,klo:khi)
            y(1:2)   = sptab(ind2,node,klo:khi)
            yp(1:2)  = sptab(ind3,node,klo:khi)
            sig(1:2) = sptab(ind3+2,node,klo:khi)
            if (iWhat <= 2) then
               ! WC (1) or K (2)
               ye = my_HVAL (xe_local,x,y,yp,sig,IER); if (IER /= 0) write (124,*) IER
               if (do_ln_trans .and. iWhat == 2) ye = dexp(ye)      ! back transformation for K
            else
               ! C = dWC/dh (3) or dK/dh (4)
               dyedxe = my_HPVAL (xe_local,x,y,yp,sig,IER); if (IER /= 0) write (125,*) IER
               if (do_ln_trans .and. iWhat == 3) dyedxe =  dyedxe/(-xe+1.0d0)           ! back transformation for C
               if (do_ln_trans .and. iWhat == 4) then
                  ye = my_HVAL (xe_local,x,y,yp,sig,IER); if (IER /= 0) write (124,*) IER
                  ye = dexp(ye)                          ! back transformation for K
                  dyedxe =  dyedxe*ye/(-xe+1.0d0)        ! back transformation for dK/dh
               end if
            end if
            return
         else
            f1 = sptab(ind2,node,klo)
            f2 = sptab(ind2,node,khi)
            d1 = sptab(ind3,node,klo)
            d2 = sptab(ind3,node,khi)
         
!           coefficients to be used
            h     = x2 - x1
            delta = (f2 - f1)/h
            del1  = (d1 - delta)/h
            del2  = (d2 - delta)/h
            c2    = -(del1+del1 + del2)
            c2t2  = c2 + c2
            c3    = (del1 + del2)/h
            c3t3  = c3+c3+c3
!           distance to reference of interval
            xx = xe_local - x1
!           compute
            ye = f1 + xx*(d1 + xx*(c2 + xx*c3))
            if (do_ln_trans .and. ind2 == 3) ye = dexp(ye)      ! back transformation for K
         end if
            
      else if(inverse.eq.1)then
!     bi-section search method for finding the array entry
         klo=1 
         khi=n
2000     if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(sptab(ind2,node,k).gt.ye)then
               khi=k
            else
               klo=k
            endif
            goto 2000
         endif 
         x1 = sptab(ind1,node,klo)
         x2 = sptab(ind1,node,khi)
         f1 = sptab(ind2,node,klo)
         f2 = sptab(ind2,node,khi)
         d1 = sptab(ind3,node,klo)
         d2 = sptab(ind3,node,khi)
         h     = x2 - x1
         delta = (f2 - f1)/h
         del1  = (d1 - delta)/h
         del2  = (d2 - delta)/h
         c2    = -(del1+del1 + del2)
         c2t2  = c2 + c2
         c3    = (del1 + del2)/h
         c3t3  = c3+c3+c3

!     bi-section search method for finding xx-value
         xlow   = x1
         xhigh  = x2
         maxtry = 20
         ntry   = 1
         xtry = 0.5d0*(xlow+xhigh)
         xx = xtry - x1
         ytry = f1 + xx*(d1 + xx*(c2 + xx*c3))
         do while (dabs(ye-ytry).gt.1.0d-05 .and.ntry.lt.maxtry)
            ntry = ntry + 1
            if(ytry.lt.ye)then
               xlow = xtry
            else
               xhigh = xtry
            end if
            xtry = 0.5d0*(xlow+xhigh)
            xx = xtry - x1
            ytry = f1 + xx*(d1 + xx*(c2 + xx*c3))
         end do
         xe = xx + x1
         if (do_ln_trans) xe = -(dexp(-xe)-1.0d0)
      end if

      dyedxe = d1 + xx*(c2t2 + xx*c3t3)
      if (do_ln_trans .and. ind2 == 2) dyedxe =  dyedxe/(-xe+1.0d0)           ! back transformation for C
      if (do_ln_trans .and. ind2 == 3) dyedxe =  dyedxe*ye/(-xe+1.0d0)        ! back transformation for K
      
      return
      end

      subroutine PreProcTabulatedFunction(flag,n,x,y,dydx,sigma)
      use doTSPACK
      use TSPACK
      implicit none
      include 'arrays.fi'

      integer  i,n, flag
      real(8)  y(matab), dydx(matab), sigma(matab)
      integer  ic(2), INCFD, NWK, IERR
      real(8)  vc(2), switch, x(matab)
      real(8)  f(1,matab),d(1,matab),wk(2*matab)
      
      if (use_TSPACK) then
         B = 0.0d0; WK2 = 0.0d0; SIGMA = 0.0d0; ICFLG = 0

         per   = .false.   ! no periodic function
         NCD   = 2         ! Number of continuous derivatives at the knots (1 or 2).
         BMAX  = 1.0d30    ! user-defined value for 'infinity'
!           B(,): lower/upper bound contraints for H and HP for each interval in X; a NULL-constraint applies if B > |BMAX|
!           B(5,): defines the sign of HPP (+ = convex,; - = concave; 0 = no restriction)
         B(1,1:N) =  1.1d0*BMAX     ! upper bound on H
         B(2,1:N) = -1.1d0*BMAX     ! lower bound on H
         B(3,1:N) =  1.1d0*BMAX     ! upper bound on HP
         B(4,1:N) = -1.1d0*BMAX     ! lower bound on HP
         B(5,1:N) =  1.0d0          ! convex (+) or concave (-)

!        WC(h) and dWC/dh as YP: x = h, y = WC, dydx = dWC/dh = YP
         if (flag == 1) then
            IENDC = 1         ! End condition indicator for NCD = 2 and PER = .FALSE.; 1 = first derivative; 2 = second derivative
            dydx(1) = 0.0d0
            dydx(N) = 1.0d-7  ! 0.0d0
            call TSPBI (N, x, y, NCD, IENDC, PER, B, BMAX, LWK, WK2, dydx, sigma, ICFLG, IER)

!        K(h) and dK/dh as YP_K: x = h, y = K, dydx = dK/dh = YP_K
         else if (flag == 2) then
            IENDC_K = 2         ! End condition indicator for NCD = 2 and PER = .FALSE.; 1 = first derivative; 2 = second derivative
            dydx(1) = 0.0d0
            dydx(N) = 0.0d0
            WK2 = 0.0d0; ICFLG = 0
            call TSPBI (N, x, y, NCD, IENDC_K, PER, B, BMAX, LWK, WK2, dydx, sigma, ICFLG, IER)
         end if

      else
!---- initialize the dydx computations
!---- second derivative at x(1) is assumed to be zero
         ic(1) = 2    
         vc(1) = 0.0d0

!---- first derivative at x(n) is assumed to be 1.0d-7 (saturation)
         if(flag.eq.1)then
            ic(2) = 1
            vc(2) = 1.0d-07
         else if(flag.eq.2)then
            ic(2) = 2
            vc(2) = 0.0
         end if

         SWITCH = 0.0d0
         INCFD = 1
         NWK = 2*(n-1)
         do i=1,n
            f(1,i) = y(i)
         end do

!---- perform the dydx computations

         call DPCHIC (IC, VC, SWITCH,n,x,f,d,INCFD,WK, NWK, IERR)

!---- store results in dydx array

         do i=1,n
            dydx(i) = d(1,i)
         end do

      end if
      return
      end

      SUBROUTINE DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK,    &
     &   IERR)
!C***BEGIN PROLOGUE  DPCHIC
!C***PURPOSE  Set derivatives needed to determine a piecewise monotone
!C            piecewise cubic Hermite interpolant to given data.
!C            User control is available over boundary conditions and/or
!C            treatment of points where monotonicity switches direction.
!C***LIBRARY   SLATEC (PCHIP)
!C***CATEGORY  E1A
!C***TYPE      DOUBLE PRECISION (PCHIC-S, DPCHIC-D)
!C***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!C             PCHIP, PIECEWISE CUBIC INTERPOLATION,
!C             SHAPE-PRESERVING INTERPOLATION
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C             Lawrence Livermore National Laboratory
!C             P.O. Box 808  (L-316)
!C             Livermore, CA  94550
!C             FTS 532-4275, (510) 422-4275
!C***DESCRIPTION
!C
!C         DPCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
!C
!C     Sets derivatives needed to determine a piecewise monotone piece-
!C     wise cubic interpolant to the data given in X and F satisfying the
!C     boundary conditions specified by IC and VC.
!C
!C     The treatment of points where monotonicity switches direction is
!C     controlled by argument SWITCH.
!C
!C     To facilitate two-dimensional applications, includes an increment
!C     between successive values of the F- and D-arrays.
!C
!C     The resulting piecewise cubic Hermite function may be evaluated
!C     by DPCHFE or DPCHFD.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  IC(2), N, NWK, IERR
!C        DOUBLE PRECISION  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),
!C                          WK(NWK)
!C
!C        CALL DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!C
!C   Parameters:
!C
!C     IC -- (input) integer array of length 2 specifying desired
!C           boundary conditions:
!C           IC(1) = IBEG, desired condition at beginning of data.
!C           IC(2) = IEND, desired condition at end of data.
!C
!C           IBEG = 0  for the default boundary condition (the same as
!C                     used by DPCHIM).
!C           If IBEG.NE.0, then its sign indicates whether the boundary
!C                     derivative is to be adjusted, if necessary, to be
!C                     compatible with monotonicity:
!C              IBEG.GT.0  if no adjustment is to be performed.
!C              IBEG.LT.0  if the derivative is to be adjusted for
!C                     monotonicity.
!C
!C           Allowable values for the magnitude of IBEG are:
!C           IBEG = 1  if first derivative at X(1) is given in VC(1).
!C           IBEG = 2  if second derivative at X(1) is given in VC(1).
!C           IBEG = 3  to use the 3-point difference formula for D(1).
!C                     (Reverts to the default b.c. if N.LT.3 .)
!C           IBEG = 4  to use the 4-point difference formula for D(1).
!C                     (Reverts to the default b.c. if N.LT.4 .)
!C           IBEG = 5  to set D(1) so that the second derivative is con-
!C              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
!C              This option is somewhat analogous to the "not a knot"
!C              boundary condition provided by DPCHSP.
!C
!C          NOTES (IBEG):
!C           1. An error return is taken if ABS(IBEG).GT.5 .
!C           2. Only in case  IBEG.LE.0  is it guaranteed that the
!C              interpolant will be monotonic in the first interval.
!C              If the returned value of D(1) lies between zero and
!C              3*SLOPE(1), the interpolant will be monotonic.  This
!C              is **NOT** checked if IBEG.GT.0 .
!C           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
!C              tonicity, a warning error is returned.
!C
!C           IEND may take on the same values as IBEG, but applied to
!C           derivative at X(N).  In case IEND = 1 or 2, the value is
!C           given in VC(2).
!C
!C          NOTES (IEND):
!C           1. An error return is taken if ABS(IEND).GT.5 .
!C           2. Only in case  IEND.LE.0  is it guaranteed that the
!C              interpolant will be monotonic in the last interval.
!C              If the returned value of D(1+(N-1)*INCFD) lies between
!C              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!C              This is **NOT** checked if IEND.GT.0 .
!C           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
!C              achieve monotonicity, a warning error is returned.
!C
!C     VC -- (input) real(8)array of length 2 specifying desired boundary
!C           values, as indicated above.
!C           VC(1) need be set only if IC(1) = 1 or 2 .
!C           VC(2) need be set only if IC(2) = 1 or 2 .
!C
!C     SWITCH -- (input) indicates desired treatment of points where
!C           direction of monotonicity switches:
!C           Set SWITCH to zero if interpolant is required to be mono-
!C           tonic in each interval, regardless of monotonicity of data.
!C             NOTES:
!C              1. This will cause D to be set to zero at all switch
!C                 points, thus forcing extrema there.
!C              2. The result of using this option with the default boun-
!C                 dary conditions will be identical to using DPCHIM, but
!C                 will generally cost more compute time.
!C                 This option is provided only to facilitate comparison
!C                 of different switch and/or boundary conditions.
!C           Set SWITCH nonzero to use a formula based on the 3-point
!C              difference formula in the vicinity of switch points.
!C           If SWITCH is positive, the interpolant on each interval
!C              containing an extremum is controlled to not deviate from
!C              the data by more than SWITCH*DFLOC, where DFLOC is the
!C              maximum of the change of F on this interval and its two
!C              immediate neighbors.
!C           If SWITCH is negative, no such control is to be imposed.
!C
!C     N -- (input) number of data points.  (Error return if N.LT.2 .)
!C
!C     X -- (input) real(8)array of independent variable values.  The
!C           elements of X must be strictly increasing:
!C                X(I-1) .LT. X(I),  I = 2(1)N.
!C           (Error return if not.)
!C
!C     F -- (input) real(8)array of dependent variable values to be
!C           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!C           X(I).
!C
!C     D -- (output) real(8)array of derivative values at the data
!C           points.  These values will determine a monotone cubic
!C           Hermite function on each subinterval on which the data
!C           are monotonic, except possibly adjacent to switches in
!C           monotonicity. The value corresponding to X(I) is stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C           No other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in F and D.
!C           This argument is provided primarily for 2-D applications.
!C           (Error return if  INCFD.LT.1 .)
!C
!C     WK -- (scratch) real(8)array of working storage.  The user may
!C           wish to know that the returned values are:
!C              WK(I)     = H(I)     = X(I+1) - X(I) ;
!C              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!C           for  I = 1(1)N-1.
!C
!C     NWK -- (input) length of work array.
!C           (Error return if  NWK.LT.2*(N-1) .)
!C
!C     IERR -- (output) error flag.
!C           Normal return:
!C              IERR = 0  (no errors).
!C           Warning errors:
!C              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!C                        monotonicity.
!C              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!C                        adjusted for monotonicity.
!C              IERR = 3  if both of the above are true.
!C           "Recoverable" errors:
!C              IERR = -1  if N.LT.2 .
!C              IERR = -2  if INCFD.LT.1 .
!C              IERR = -3  if the X-array is not strictly increasing.
!C              IERR = -4  if ABS(IBEG).GT.5 .
!C              IERR = -5  if ABS(IEND).GT.5 .
!C              IERR = -6  if both of the above are true.
!C              IERR = -7  if NWK.LT.2*(N-1) .
!C             (The D-array has not been changed in any of these cases.)
!C               NOTE:  The above errors are checked in the order listed,
!C                   and following arguments have **NOT** been validated.
!C
!C***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
!C                 Package, Report UCRL-87285, Lawrence Livermore Natio-
!C                 nal Laboratory, July 1982.  [Poster presented at the
!C                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
!C               2. F. N. Fritsch and J. Butland, A method for construc-
!C                 ting local monotone piecewise cubic interpolants, SIAM
!C                 Journal on Scientific and Statistical Computing 5, 2
!C                 (June 1984), pp. 300-304.
!C               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!C                 cubic interpolation, SIAM Journal on Numerical Ana-
!C                 lysis 17, 2 (April 1980), pp. 238-246.
!!
!C***ROUTINES CALLED  DPCHCE, DPCHCI, DPCHCS, XERMSG
!!
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820804  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   870813  Updated Reference 2.
!C   890206  Corrected XERROR calls.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890703  Corrected category record.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   920429  Revised format and order of references.  (WRB,FNF)
!C***END PROLOGUE  DPCHIC
!C  Programming notes:
!C
!C     To produce a single precision version, simply:
!C        a. Change DPCHIC to PCHIC wherever it occurs,
!C        b. Change DPCHCE to PCHCE wherever it occurs,
!C        c. Change DPCHCI to PCHCI wherever it occurs,
!C        d. Change DPCHCS to PCHCS wherever it occurs,
!C        e. Change the double precision declarations to real, and
!C        f. Change the constant  ZERO  to single precision.
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  IC(2), N, INCFD, NWK, IERR
      DOUBLE PRECISION  VC(2), SWITCH, X(*), F(INCFD,*), D(INCFD,*),    &
     & WK(NWK)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, IBEG, IEND, NLESS1
      DOUBLE PRECISION  vsmall
      DATA  vsmall/1.d-15/
!C
!C  VALIDITY-CHECK ARGUMENTS.
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHIC
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF (ABS(IBEG) .GT. 5)  IERR = IERR - 1
      IF (ABS(IEND) .GT. 5)  IERR = IERR - 2
      IF (IERR .LT. 0)  GO TO 5004
!C
!C  FUNCTION DEFINITION IS OK -- GO ON.
!C
      NLESS1 = N - 1
      IF ( NWK .LT. 2*NLESS1 )  GO TO 5007
!C
!C  SET UP H AND SLOPE ARRAYS.
!C
      DO 20  I = 1, NLESS1
         WK(I) = X(I+1) - X(I)
         WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
   20 CONTINUE
!C
!C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!C
      IF (NLESS1 .GT. 1)  GO TO 1000
      D(1,1) = WK(2)
      D(1,N) = WK(2)
      GO TO 3000
!C
!C  NORMAL CASE  (N .GE. 3) .
!C
 1000 CONTINUE
!C
!C  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!C
!C     --------------------------------------
      CALL DPCHCI (N, WK(1:N-1), WK(N:2*N-2), D, INCFD)  
!C     --------------------------------------
!C
!C  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!C
      IF (dabs(SWITCH).lt.vsmall)  GO TO 3000

!C     ----------------------------------------------------
      CALL DPCHCS (SWITCH, N, WK(1:N-1), WK(N:2*N-2), D, INCFD, IERR)
!C     ----------------------------------------------------
      IF (IERR .NE. 0)  GO TO 5008
!C
!C  SET END CONDITIONS.
!C
 3000 CONTINUE
      IF ( (IBEG.EQ.0) .AND. (IEND.EQ.0) )  GO TO 5000
!C     -------------------------------------------------------
      CALL DPCHCE (IC, VC, N, X, WK(1:N-1), WK(N:2*N-2), D, INCFD, IERR)
!C     -------------------------------------------------------
      IF (IERR .LT. 0)  GO TO 5009
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      RETURN
!C
!C  ERROR RETURNS.
!C
 5001 CONTINUE
!C     N.LT.2 RETURN.
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHIC',
!     +   'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
!C
 5002 CONTINUE
!C     INCFD.LT.1 RETURN.
      IERR = -2
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'INCREMENT LESS THAN ONE', IERR,
!     +   1)
      RETURN
!C
 5003 CONTINUE
!C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
!      CALL XERMSG ('SLATEC', 'DPCHIC',
!     +   'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
      RETURN
!C
 5004 CONTINUE
!C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'IC OUT OF RANGE', IERR, 1)
      RETURN
!C
 5007 CONTINUE
!C     NWK .LT. 2*(N-1)  RETURN.
      IERR = -7
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'WORK ARRAY TOO SMALL', IERR, 1)
      RETURN
!C
 5008 CONTINUE
!C     ERROR RETURN FROM DPCHCS.
      IERR = -8
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCS',
!     +   IERR, 1)
      RETURN
!C
 5009 CONTINUE
!C     ERROR RETURN FROM DPCHCE.
!C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
!      CALL XERMSG ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCE',
!     +   IERR, 1)
      RETURN
!C------------- LAST LINE OF DPCHIC FOLLOWS -----------------------------
      END

      SUBROUTINE DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!C***BEGIN PROLOGUE  DPCHCE
!C***SUBSIDIARY
!C***PURPOSE  Set boundary conditions for DPCHIC
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHCE-S, DPCHCE-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C          DPCHCE:  DPCHIC End Derivative Setter.
!C
!C    Called by DPCHIC to set end derivatives as requested by the user.
!C    It must be called after interior derivative values have been set.
!C                      -----
!C
!C    To facilitate two-dimensional applications, includes an increment
!C    between successive values of the D-array.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  IC(2), N, IERR
!C        DOUBLE PRECISION  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!C
!C        CALL  DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!C
!C   Parameters:
!C
!C     IC -- (input) integer array of length 2 specifying desired
!C           boundary conditions:
!C           IC(1) = IBEG, desired condition at beginning of data.
!C           IC(2) = IEND, desired condition at end of data.
!C           ( see prologue to DPCHIC for details. )
!C
!C     VC -- (input) real(8)array of length 2 specifying desired boundary
!C           values.  VC(1) need be set only if IC(1) = 2 or 3 .
!C                    VC(2) need be set only if IC(2) = 2 or 3 .
!C
!C     N -- (input) number of data points.  (assumes N.GE.2)
!C
!C     X -- (input) real(8)array of independent variable values.  (the
!C           elements of X are assumed to be strictly increasing.)
!C
!C     H -- (input) real(8)array of interval lengths.
!C     SLOPE -- (input) real(8)array of data slopes.
!C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!C                  H(I) =  X(I+1)-X(I),
!C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!C
!C     D -- (input) real(8)array of derivative values at the data points.
!C           The value corresponding to X(I) must be stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C          (output) the value of D at X(1) and/or X(N) is changed, if
!C           necessary, to produce the requested boundary conditions.
!C           no other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in D.
!C           This argument is provided primarily for 2-D applications.
!C
!C     IERR -- (output) error flag.
!C           Normal return:
!C              IERR = 0  (no errors).
!C           Warning errors:
!C              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!C                        monotonicity.
!C              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!C                        adjusted for monotonicity.
!C              IERR = 3  if both of the above are true.
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS.
!C
!C***SEE ALSO  DPCHIC
!C***ROUTINES CALLED  DPCHDF, DPCHST, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   890206  Corrected XERROR calls.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR section in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHCE
!C
!C  Programming notes:
!C     1. The function DPCHST(ARG1,ARG2)  is assumed to return zero if
!C        either argument is zero, +1 if they are of the same sign, and
!C        -1 if they are of opposite sign.
!C     2. One could reduce the number of arguments and amount of local
!C        storage, at the expense of reduced code clarity, by passing in
!C        the array WK (rather than splitting it into H and SLOPE) and
!C        increasing its length enough to incorporate STEMP and XTEMP.
!C     3. The two monotonicity checks only use the sufficient conditions.
!C        Thus, it is possible (but unlikely) for a boundary condition to
!C        be changed, even though the original interpolant was monotonic.
!C        (At least the result is a continuous function of the data.)
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  IC(2), N, INCFD, IERR
      DOUBLE PRECISION  VC(2), X(*), H(*), SLOPE(*), D(INCFD,*)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  IBEG, IEND, IERF, INDEX, J, K
      DOUBLE PRECISION  HALF, STEMP(3), THREE, TWO, XTEMP(4), ZERO
      SAVE ZERO, HALF, TWO, THREE
      DOUBLE PRECISION  DPCHDF, DPCHST, vsmall
!C
!C  INITIALIZE.
!C
      DATA  ZERO /0.D0/,  HALF/.5D0/,  TWO/2.D0/, THREE/3.D0/
      DATA  vsmall /1.D-15/
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHCE
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
!C
!C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!C
      IF ( ABS(IBEG).GT.N )  IBEG = 0
      IF ( ABS(IEND).GT.N )  IEND = 0
!C
!C  TREAT BEGINNING BOUNDARY CONDITION.
!C
      IF (IBEG .EQ. 0)  GO TO 2000
      K = ABS(IBEG)
      IF (K .EQ. 1)  THEN
!C        BOUNDARY VALUE PROVIDED.
         D(1,1) = VC(1)
      ELSE IF (K .EQ. 2)  THEN
!C        BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,1) = HALF*( (THREE*SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
      ELSE IF (K .LT. 5)  THEN
!C        USE K-POINT DERIVATIVE FORMULA.
!C        PICK UP FIRST K POINTS, IN REVERSE ORDER.
         DO 10  J = 1, K
            INDEX = K-J+1
!C           INDEX RUNS FROM K DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1)
   10    CONTINUE
!C                 -----------------------------
         D(1,1) = DPCHDF (K, XTEMP, STEMP, IERF)
!C                 -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
!C        USE 'NOT A KNOT' CONDITION.
         D(1,1) = ( THREE*(H(1)*SLOPE(2) + H(2)*SLOPE(1))               &
     &             - TWO*(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
      ENDIF
!C
      IF (IBEG .GT. 0)  GO TO 2000
!C
!C  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!C
      IF (dabs(SLOPE(1)).lt.vsmall)  THEN
         IF (dabs(D(1,1)).gt.ZERO)  THEN
            D(1,1) = ZERO
            IERR = IERR + 1
         ENDIF
      ELSE IF ( DPCHST(D(1,1),SLOPE(1)) .LT. ZERO)  THEN
         D(1,1) = ZERO
         IERR = IERR + 1
      ELSE IF ( DABS(D(1,1)) .GT. THREE*DABS(SLOPE(1)) )  THEN
         D(1,1) = THREE*SLOPE(1)
         IERR = IERR + 1
      ENDIF
!C
!C  TREAT END BOUNDARY CONDITION.
!C
 2000 CONTINUE
      IF (IEND .EQ. 0)  GO TO 5000
      K = ABS(IEND)
      IF (K .EQ. 1)  THEN
!C        BOUNDARY VALUE PROVIDED.
         D(1,N) = VC(2)
      ELSE IF (K .EQ. 2)  THEN
!C        BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,N) = HALF*( (THREE*SLOPE(N-1) - D(1,N-1)) +                &
     &                                           HALF*VC(2)*H(N-1) )
      ELSE IF (K .LT. 5)  THEN
!C        USE K-POINT DERIVATIVE FORMULA.
!C        PICK UP LAST K POINTS.
         DO 2010  J = 1, K
            INDEX = N-K+J
!C           INDEX RUNS FROM N+1-K UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX)
 2010    CONTINUE
!C                 -----------------------------
         D(1,N) = DPCHDF (K, XTEMP, STEMP, IERF)
!C                 -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
!C        USE 'NOT A KNOT' CONDITION.
         D(1,N) = ( THREE*(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1))       &
     &             - TWO*(H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) )   &
     &                                                         / H(N-2)
      ENDIF
!C
      IF (IEND .GT. 0)  GO TO 5000
!C
!C  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!C

      IF (dabs(SLOPE(N-1)).lt.vsmall)  THEN
         IF (dabs(D(1,N)).gt.ZERO)  THEN
            D(1,N) = ZERO
            IERR = IERR + 2
         ENDIF
      ELSE IF ( DPCHST(D(1,N),SLOPE(N-1)) .LT. ZERO)  THEN
         D(1,N) = ZERO
         IERR = IERR + 2
      ELSE IF ( DABS(D(1,N)) .GT. THREE*DABS(SLOPE(N-1)) )  THEN
         D(1,N) = THREE*SLOPE(N-1)
         IERR = IERR + 2
      ENDIF
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      RETURN
!C
!C  ERROR RETURN.
!C
 5001 CONTINUE
!C     ERROR RETURN FROM DPCHDF.
!C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHCE', 'ERROR RETURN FROM DPCHDF',
!     +   IERR, 1)
      RETURN
!C------------- LAST LINE OF DPCHCE FOLLOWS -----------------------------
      END

      SUBROUTINE DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!C***BEGIN PROLOGUE  DPCHCS
!C***SUBSIDIARY
!C***PURPOSE  Adjusts derivative values for DPCHIC
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHCS-S, DPCHCS-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C         DPCHCS:  DPCHIC Monotonicity Switch Derivative Setter.
!C
!C     Called by  DPCHIC  to adjust the values of D in the vicinity of a
!C     switch in direction of monotonicity, to produce a more "visually
!C     pleasing" curve than that given by  DPCHIM .
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  N, IERR
!C        DOUBLE PRECISION  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!C
!C        CALL  DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!C
!C   Parameters:
!C
!C     SWITCH -- (input) indicates the amount of control desired over
!C           local excursions from data.
!C
!C     N -- (input) number of data points.  (assumes N.GT.2 .)
!C
!C     H -- (input) real(8)array of interval lengths.
!C     SLOPE -- (input) real(8)array of data slopes.
!C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!C                  H(I) =  X(I+1)-X(I),
!C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!C
!C     D -- (input) real(8)array of derivative values at the data points,
!C           as determined by DPCHCI.
!C          (output) derivatives in the vicinity of switches in direction
!C           of monotonicity may be adjusted to produce a more "visually
!C           pleasing" curve.
!C           The value corresponding to X(I) is stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C           No other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in D.
!C           This argument is provided primarily for 2-D applications.
!C
!C     IERR -- (output) error flag.  should be zero.
!C           If negative, trouble in DPCHSW.  (should never happen.)
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS, MAX, MIN.
!C
!C***SEE ALSO  DPCHIC
!C***ROUTINES CALLED  DPCHST, DPCHSW
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820617  Redesigned to (1) fix  problem with lack of continuity
!C           approaching a flat-topped peak (2) be cleaner and
!C           easier to verify.
!C           Eliminated subroutines PCHSA and PCHSX in the process.
!C   820622  1. Limited fact to not exceed one, so computed D is a
!C             convex combination of DPCHCI value and DPCHSD value.
!C           2. Changed fudge from 1 to 4 (based on experiments).
!C   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected conversion to double precision.
!C   870813  Minor cosmetic changes.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891006  Modified spacing in computation of DFLOC.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR section in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHCS
!C
!C  Programming notes:
!C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!C        either argument is zero, +1 if they are of the same sign, and
!C        -1 if they are of opposite sign.
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  N, INCFD, IERR
      DOUBLE PRECISION  SWITCH, H(*), SLOPE(*), D(INCFD,*)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, INDX, K, NLESS1
      DOUBLE PRECISION  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE,    &
     &      SLMAX, WTAVE(2), ZERO
      SAVE ZERO, ONE, FUDGE
      DOUBLE PRECISION  DPCHST
!C
!C  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
!C
!      DOUBLE PRECISION  DPCHSD , S1, S2, H1, H2
!      DPCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
!C
!C  INITIALIZE.
!C
      DATA  ZERO /0.D0/,  ONE/1.D0/
      DATA  FUDGE /4.D0/
!C***FIRST EXECUTABLE STATEMENT  DPCHCS
      IERR = 0
      NLESS1 = N - 1
!C
!C  LOOP OVER SEGMENTS.
!C
      DO 900  I = 2, NLESS1
! obsolete arithmetic if         IF ( DPCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900
         if (DPCHST(SLOPE(I-1),SLOPE(I)) < 0.d0) then
            goto 100
         else if (DPCHST(SLOPE(I-1),SLOPE(I)) > 0.d0) then
            goto 900
         else
            goto 300
         end if
!C             --------------------------
!C
  100    CONTINUE
!C
!C....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
!C
!C           DO NOT CHANGE D IF 'UP-DOWN-UP'.
            IF (I .GT. 2)  THEN
               IF ( DPCHST(SLOPE(I-2),SLOPE(I)) .GT. ZERO)  GO TO 900
!C                   --------------------------
            ENDIF
            IF (I .LT. NLESS1)  THEN
               IF ( DPCHST(SLOPE(I+1),SLOPE(I-1)) .GT. ZERO)  GO TO 900
!C                   ----------------------------
            ENDIF
!C
!C   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
!C
            DEXT = DPCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
!C
!C   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!C
! obsolete arithmetic if            IF ( DPCHST(DEXT, SLOPE(I-1)) )  200, 900, 250
            if (DPCHST(DEXT, SLOPE(I-1)) < 0.d0) then
               goto 200
            else if (DPCHST(DEXT, SLOPE(I-1)) > 0.d0) then
               goto 250
            else
               goto 900
            end if
!C                -----------------------
!C
  200       CONTINUE
!C              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
!C                        EXTREMUM IS IN (X(I-1),X(I)).
               K = I-1
!C              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
               WTAVE(2) = DEXT
               IF (K .GT. 1)                                            &
     &            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
               GO TO 400
!C
  250       CONTINUE
!C              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
!C                        EXTREMUM IS IN (X(I),X(I+1)).
               K = I
!C              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
               WTAVE(1) = DEXT
               IF (K .LT. NLESS1)                                       &
     &            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
               GO TO 400
!C
  300    CONTINUE
!C
!C....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
!C                     CHECK FOR FLAT-TOPPED PEAK .......................
!C
            IF (I .EQ. NLESS1)  GO TO 900
            IF ( DPCHST(SLOPE(I-1), SLOPE(I+1)) .GE. ZERO)  GO TO 900
!C                -----------------------------
!C
!C           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
            K = I
!C           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
!C
  400    CONTINUE
!C
!C....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!C        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!C           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
!C                    IF K.GT.1
!C           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
!C                    IF K.LT.N-1
!C
         SLMAX = DABS(SLOPE(K))
         IF (K .GT. 1)    SLMAX = MAX( SLMAX, DABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) SLMAX = MAX( SLMAX, DABS(SLOPE(K+1)) )
!C
         IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX
         DEL(2) = SLOPE(K) / SLMAX
         IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX
!C
         IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN
!C           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
            FACT = FUDGE* DABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
            D(1,K) = D(1,K) + MIN(FACT,ONE)*(WTAVE(1) - D(1,K))
            FACT = FUDGE* DABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
            D(1,K+1) = D(1,K+1) + MIN(FACT,ONE)*(WTAVE(2) - D(1,K+1))
         ELSE
!C           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
!C                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
            FACT = FUDGE* DABS(DEL(2))
            D(1,I) = MIN(FACT,ONE) * WTAVE(I-K+1)
!C              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
!C                        I-K+1 = 2 IF K=I-1(=1).
         ENDIF
!C
!C
!C....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!C
         IF (SWITCH .LE. ZERO)  GO TO 900
!C
         DFLOC = H(K)*DABS(SLOPE(K))
         IF (K .GT. 1)    DFLOC = MAX( DFLOC, H(K-1)*DABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) DFLOC = MAX( DFLOC, H(K+1)*DABS(SLOPE(K+1)) )
         DFMX = SWITCH*DFLOC
         INDX = I-K+1
!C        INDX = 1 IF K=I, 2 IF K=I-1.
!C        ---------------------------------------------------------------
         CALL DPCHSW(DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)
!C        ---------------------------------------------------------------
         IF (IERR .NE. 0)  RETURN
!C
!C....... END OF SEGMENT LOOP.
!C
  900 CONTINUE
!C
      RETURN
!C------------- LAST LINE OF DPCHCS FOLLOWS -----------------------------
   contains
   
   double precision function DPCHSD (S1,S2,H1,H2)
   double precision S1,S2,H1,H2
   DPCHSD = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
   return
   end function DPCHSD
         
      END

      SUBROUTINE DPCHCI (N, H, SLOPE, D, INCFD)
!C***BEGIN PROLOGUE  DPCHCI
!C***SUBSIDIARY
!C***PURPOSE  Set interior derivatives for DPCHIC
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHCI-S, DPCHCI-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C          DPCHCI:  DPCHIC Initial Derivative Setter.
!C
!C    Called by DPCHIC to set derivatives needed to determine a monotone
!C    piecewise cubic Hermite interpolant to the data.
!C
!C    Default boundary conditions are provided which are compatible
!C    with monotonicity.  If the data are only piecewise monotonic, the
!C    interpolant will have an extremum at each point where monotonicity
!C    switches direction.
!C
!C    To facilitate two-dimensional applications, includes an increment
!C    between successive values of the D-array.
!C
!C    The resulting piecewise cubic Hermite function should be identical
!C    (within roundoff error) to that produced by DPCHIM.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        PARAMETER  (INCFD = ...)
!C        INTEGER  N
!C        DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
!C
!C        CALL  DPCHCI (N, H, SLOPE, D, INCFD)
!C
!C   Parameters:
!C
!C     N -- (input) number of data points.
!C           If N=2, simply does linear interpolation.
!C
!C     H -- (input) real(8)array of interval lengths.
!C     SLOPE -- (input) real(8)array of data slopes.
!C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!C                  H(I) =  X(I+1)-X(I),
!C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!C
!C     D -- (output) real(8)array of derivative values at data points.
!C           If the data are monotonic, these values will determine a
!C           a monotone cubic Hermite function.
!C           The value corresponding to X(I) is stored in
!C                D(1+(I-1)*INCFD),  I=1(1)N.
!C           No other entries in D are changed.
!C
!C     INCFD -- (input) increment between successive values in D.
!C           This argument is provided primarily for 2-D applications.
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS, MAX, MIN.
!C
!C***SEE ALSO  DPCHIC
!C***ROUTINES CALLED  DPCHST
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820601  Modified end conditions to be continuous functions of
!C           data when monotonicity switches in next interval.
!C   820602  1. Modified formulas so end conditions are less prone
!C             to over/underflow problems.
!C           2. Minor modification to HSUM calculation.
!C   820805  Converted to SLATEC library version.
!C   870813  Minor cosmetic changes.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR section in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHCI
!C
!C  Programming notes:
!C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!C        either argument is zero, +1 if they are of the same sign, and
!C        -1 if they are of opposite sign.
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  N, INCFD
      DOUBLE PRECISION  H(*), SLOPE(*), D(INCFD,*)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, NLESS1
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM,     &
     &      HSUMT3, THREE, W1, W2, ZERO
      SAVE ZERO, THREE
      DOUBLE PRECISION  DPCHST
!C
!C  INITIALIZE.
!C
      DATA  ZERO /0.D0/, THREE/3.D0/
!C***FIRST EXECUTABLE STATEMENT  DPCHCI
      NLESS1 = N - 1
      DEL1 = SLOPE(1)
!C
!C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
!C
!C  NORMAL CASE  (N .GE. 3).
!C
   10 CONTINUE
      DEL2 = SLOPE(2)
!C
!C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!C     SHAPE-PRESERVING.
!C
      HSUM = H(1) + H(2)
      W1 = (H(1) + HSUM)/HSUM
      W2 = -H(1)/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (DABS(D(1,1)) .GT. DABS(DMAX))  D(1,1) = DMAX
      ENDIF
!C
!C  LOOP THROUGH INTERIOR POINTS.
!C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
!C
         HSUM = H(I-1) + H(I)
         DEL1 = DEL2
         DEL2 = SLOPE(I)
   40    CONTINUE
!C
!C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!C
         D(1,I) = ZERO
         IF ( DPCHST(DEL1,DEL2) .LE. ZERO)  GO TO 50
!C
!C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!C
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H(I-1))/HSUMT3
         W2 = (HSUM + H(I)  )/HSUMT3
         DMAX = MAX( DABS(DEL1), DABS(DEL2) )
         DMIN = MIN( DABS(DEL1), DABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!C
   50 CONTINUE
!C
!C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!C     SHAPE-PRESERVING.
!C
      W1 = -H(N-1)/HSUM
      W2 = (H(N-1) + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (DABS(D(1,N)) .GT. DABS(DMAX))  D(1,N) = DMAX
      ENDIF
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      RETURN
!C------------- LAST LINE OF DPCHCI FOLLOWS -----------------------------
      END

      DOUBLE PRECISION FUNCTION DPCHDF (K, X, S, IERR)
!C***BEGIN PROLOGUE  DPCHDF
!C***SUBSIDIARY
!C***PURPOSE  Computes divided differences for DPCHCE and DPCHSP
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHDF-S, DPCHDF-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C          DPCHDF:   DPCHIP Finite Difference Formula
!C
!C     Uses a divided difference formulation to compute a K-point approx-
!C     imation to the derivative at X(K) based on the data in X and S.
!C
!C     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
!C     derivative approximations.
!C
!C ----------------------------------------------------------------------
!C
!C     On input:
!C        K      is the order of the desired derivative approximation.
!C               K must be at least 3 (error return if not).
!C        X      contains the K values of the independent variable.
!C               X need not be ordered, but the values **MUST** be
!C               distinct.  (Not checked here.)
!C        S      contains the associated slope values:
!C                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!C               (Note that S need only be of length K-1.)
!C
!C     On return:
!C        S      will be destroyed.
!C        IERR   will be set to -1 if K.LT.2 .
!C        DPCHDF  will be set to the desired derivative approximation if
!C               IERR=0 or to zero if IERR=-1.
!C
!C ----------------------------------------------------------------------
!C
!C***SEE ALSO  DPCHCE, DPCHSP
!C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
!C                 Verlag, New York, 1978, pp. 10-16.
!C***ROUTINES CALLED  XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   820503  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   870813  Minor cosmetic changes.
!C   890206  Corrected XERROR calls.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890411  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!C   920429  Revised format and order of references.  (WRB,FNF)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHDF
!C
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  K, IERR
      DOUBLE PRECISION  X(K), S(K)
!C
!C  DECLARE LOCAL VARIABLES.
!C
      INTEGER  I, J
      DOUBLE PRECISION  VALUE, ZERO
      SAVE ZERO
      DATA  ZERO /0.D0/
!C
!C  CHECK FOR LEGAL VALUE OF K.
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHDF
      IF (K .LT. 3)  GO TO 5001
!C
!C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
!C
!C  EVALUATE DERIVATIVE AT X(K).
!C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
!C
!C  NORMAL RETURN.
!C
      IERR = 0
      DPCHDF = VALUE
      RETURN
!C
!C  ERROR RETURN.
!C
 5001 CONTINUE
!C     K.LT.3 RETURN.
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHDF', 'K LESS THAN THREE', IERR, 1)
      DPCHDF = ZERO
      RETURN
!C------------- LAST LINE OF DPCHDF FOLLOWS -----------------------------
      END


      DOUBLE PRECISION FUNCTION DPCHST (ARG1, ARG2)
!C***BEGIN PROLOGUE  DPCHST
!C***SUBSIDIARY
!C***PURPOSE  DPCHIP Sign-Testing Routine
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHST-S, DPCHST-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C         DPCHST:  DPCHIP Sign-Testing Routine.
!C
!C
!C     Returns:
!C        -1. if ARG1 and ARG2 are of opposite sign.
!C         0. if either argument is zero.
!C        +1. if ARG1 and ARG2 are of the same sign.
!C
!C     The object is to do this without multiplying ARG1*ARG2, to avoid
!C     possible over/underflow problems.
!C
!C  Fortran intrinsics used:  SIGN.
!C
!C***SEE ALSO  DPCHCE, DPCHCI, DPCHCS, DPCHIM
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   811103  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870813  Minor cosmetic changes.
!C   890411  Added SAVE statements (Vers. 3.2).
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHST
!C
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      DOUBLE PRECISION  ARG1, ARG2
!C
!C  DECLARE LOCAL VARIABLES.
!C
      DOUBLE PRECISION  ONE, ZERO, vsmall
      SAVE ZERO, ONE
      DATA  ZERO /0.D0/,  vsmall/1.D-15/,  ONE/1.D0/
      
!C
!C  PERFORM THE TEST.
!C
!C***FIRST EXECUTABLE STATEMENT  DPCHST
      DPCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((dabs(ARG1).lt.vsmall) .OR. (dabs(ARG2).lt.vsmall))           &
     &     DPCHST = ZERO
      
!C
      RETURN
!C------------- LAST LINE OF DPCHST FOLLOWS -----------------------------
      END


      SUBROUTINE DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!C***BEGIN PROLOGUE  DPCHSW
!C***SUBSIDIARY
!C***PURPOSE  Limits excursion from data for DPCHCS
!C***LIBRARY   SLATEC (PCHIP)
!C***TYPE      DOUBLE PRECISION (PCHSW-S, DPCHSW-D)
!C***AUTHOR  Fritsch, F. N., (LLNL)
!C***DESCRIPTION
!C
!C         DPCHSW:  DPCHCS Switch Excursion Limiter.
!C
!C     Called by  DPCHCS  to adjust D1 and D2 if necessary to insure that
!C     the extremum on this interval is not further than DFMAX from the
!C     extreme data value.
!C
!C ----------------------------------------------------------------------
!C
!C  Calling sequence:
!C
!C        INTEGER  IEXTRM, IERR
!C        DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!C
!C        CALL  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!C
!C   Parameters:
!C
!C     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!C           the cubic determined by derivative values D1,D2.  (assumes
!C           DFMAX.GT.0.)
!C
!C     IEXTRM -- (input) index of the extreme data value.  (assumes
!C           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
!C
!C     D1,D2 -- (input) derivative values at the ends of the interval.
!C           (Assumes D1*D2 .LE. 0.)
!C          (output) may be modified if necessary to meet the restriction
!C           imposed by DFMAX.
!C
!C     H -- (input) interval length.  (Assumes  H.GT.0.)
!C
!C     SLOPE -- (input) data slope on the interval.
!C
!C     IERR -- (output) error flag.  should be zero.
!C           If IERR=-1, assumption on D1 and D2 is not satisfied.
!C           If IERR=-2, quadratic equation locating extremum has
!C                       negative discriminant (should never occur).
!C
!C    -------
!C    WARNING:  This routine does no validity-checking of arguments.
!C    -------
!C
!C  Fortran intrinsics used:  ABS, SIGN, SQRT.
!C
!C  ***SEE ALSO  DPCHCS
!C***ROUTINES CALLED  D1MACH, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   820218  DATE WRITTEN
!C   820805  Converted to SLATEC library version.
!C   870707  Corrected XERROR calls for d.p. name(s).
!C   870707  Replaced DATA statement for SMALL with a use of D1MACH.
!C   870813  Minor cosmetic changes.
!C   890206  Corrected XERROR calls.
!C   890411  1. Added SAVE statements (Vers. 3.2).
!C           2. Added DOUBLE PRECISION declaration for D1MACH.
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900328  Added TYPE section.  (WRB)
!C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!C   920526  Eliminated possible divide by zero problem.  (FNF)
!C   930503  Improved purpose.  (FNF)
!C***END PROLOGUE  DPCHSW
!C
!C**End
!C
!C  DECLARE ARGUMENTS.
!C
      INTEGER  IEXTRM, IERR
      DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!C
!C  DECLARE LOCAL VARIABLES.
!C
      DOUBLE PRECISION  CP, FACT, HPHI, LAMBDA, NU, ONE, PHI, RADCAL,   &
     &                  RHO, SIGMA, SMALL, THAT, THIRD, THREE, TWO, ZERO
      SAVE ZERO, ONE, TWO, THREE, FACT
      SAVE THIRD
      DOUBLE PRECISION  vsmall
!C
      DATA  ZERO /0.D0/,  ONE /1.D0/,  TWO /2.D0/, THREE /3.D0/,        &
     &      FACT /100.D0/, vsmall/1.d-15/
!C        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
      DATA  THIRD /0.33333D0/
!C
!C  NOTATION AND GENERAL REMARKS.
!C
!C     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!C     LAMBDA IS THE RATIO OF D2 TO D1.
!C     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!C     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!C           WHERE  THAT = (XHAT - X1)/H .
!C        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!C     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
!C
!C      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
!C***FIRST EXECUTABLE STATEMENT  DPCHSW
!      SMALL = FACT*D1MACH(4)
      SMALL = FACT*1.d-16
!C
!C  DO MAIN CALCULATION.
!C
      IF (dabs(D1) .lt. vsmall)  THEN  
!C
!C        SPECIAL CASE -- D1.EQ.ZERO .
!C
!C          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
         IF (dabs(D2) .lt. vsmall)  GO TO 5001
!C
         RHO = SLOPE/D2
!C          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
         IF (RHO .GE. THIRD)  GO TO 5000
         THAT = (TWO*(THREE*RHO-ONE)) / (THREE*(TWO*RHO-ONE))
         PHI = THAT**2 * ((THREE*RHO-ONE)/THREE)
!C
!C          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!C
!C          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         HPHI = H * DABS(PHI)
         IF (HPHI*DABS(D2) .GT. DFMAX)  THEN
!C           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
            D2 = SIGN (DFMAX/HPHI, D2)
         ENDIF
      ELSE
!C
         RHO = SLOPE/D1
         LAMBDA = -D2/D1
         IF (dabs(D2).lt.vsmall)  THEN
!C
!C           SPECIAL CASE -- D2.EQ.ZERO .
!C
!C             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
            IF (RHO .GE. THIRD)  GO TO 5000
            CP = TWO - THREE*RHO
            NU = ONE - TWO*RHO
            THAT = ONE / (THREE*NU)
         ELSE
            IF (LAMBDA .LE. ZERO)  GO TO 5001
!C
!C           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!C
            NU = ONE - LAMBDA - TWO*RHO
            SIGMA = ONE - RHO
            CP = NU + SIGMA
            IF (DABS(NU) .GT. SMALL)  THEN
               RADCAL = (NU - (TWO*RHO+ONE))*NU + SIGMA**2
               IF (RADCAL .LT. ZERO)  GO TO 5002
               THAT = (CP - DSQRT(RADCAL)) / (THREE*NU)
            ELSE
               THAT = ONE/(TWO*SIGMA)
            ENDIF
         ENDIF
         PHI = THAT*((NU*THAT - CP)*THAT + ONE)
!C
!C          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!C
!C          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         HPHI = H * DABS(PHI)
         IF (HPHI*DABS(D1) .GT. DFMAX)  THEN
!C           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
            D1 = SIGN (DFMAX/HPHI, D1)
            D2 = -LAMBDA*D1
         ENDIF
      ENDIF
!C
!C  NORMAL RETURN.
!C
 5000 CONTINUE
      IERR = 0
      RETURN
!C
!C  ERROR RETURNS.
!C
 5001 CONTINUE
!C     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
      IERR = -1
!      CALL XERMSG ('SLATEC', 'DPCHSW', 'D1 AND/OR D2 INVALID', IERR, 1)
      RETURN
!C
 5002 CONTINUE
!C     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
      IERR = -2
!      CALL XERMSG ('SLATEC', 'DPCHSW', 'NEGATIVE RADICAL', IERR, 1)
      RETURN
!C------------- LAST LINE OF DPCHSW FOLLOWS -----------------------------
   END
   