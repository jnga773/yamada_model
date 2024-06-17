!==============================================================================!
!                           AUTO-07: PHASE RESETTING                           !
!==============================================================================!
! This module file contains the subroutines neeed to calculate the phase
! resetting problem of the Yamada model in AUTO.
!
! This module contains the following subroutines:
! - FUNC           : The encoded vector field of the continuation problem.
! - BCND           : Define the various boundary conditions (if necessary).
!
! The following subroutines are needed for AUTO, but are not used.
! - STPNT      : Define and set initial conditions of the vector field.
! - ICND       : Define the various integral conditions (if necessary).
!
! - FOPT       : Defines the objective function for algebraic optimisation
!                problems (if necessary).
! - PVLS       : Something for algebraic problems (which I don't need to worry
!                about).
!
! The U-vector components are as follows:
!       U(1:3)   - 'Segment 1': State space components,
!       U(4:6)   - 'Segment 1': Perpendicular components,
!       U(7:9)   - 'Segment 2': State space components,
!       U(10:12) - 'Segment 2': Perpendicular components,
!       U(13:15) - 'Segment 3': State space components,
!       U(16:18) - 'Segment 4': State space components,
!
! The parameters (PAR) are as follows:
!       PAR(1:4) - Parameters of the Yamada model (gamma, A, B, and a),
!       PAR(5)   - The period of the periodic orbit (T),
!       PAR(6)   - Integer multiplier for the period of the segments (k),
!       PAR(7)   - The stable eigenvalue of the Floquet bundle (mu_s),
!       PAR(8)   - Distance from perturbed segment to \Gamma (eta),
!       PAR(9)   - Phase on \Gamma at which the perturbation starts (theta_old),
!       PAR(10)  - Phase on \Gamma at which the perturbed orbit lands (theta_new),
!       PAR(11)  - Amplitude of perurbation (A_perturb),
!       PAR(12)  - Angle at which perturbation is applied (theta_perturb),
!       PAR(13)  - Azimuthal angle at which perturbation is applied (phi_perturb).

!==============================================================================!
!                      CONTINUATION PROBLEM VECTOR FIELD                       !
!==============================================================================!
SUBROUTINE FUNC(NDIM, U, ICP, PAR, IJAC, F_out, DFDU, DFDP)
  ! Vector field encoding of the phase resetting continuation problem.
  !
  ! Input
  ! ----------
  ! NDIM   : INTEGER
  !     Dimension of the algebraic or ODE system.
  ! U   : REAL(KIND=8)
  !     State variables of the vector field.
  ! ICP    : INTEGER, ARRAY
  !     Array indicating the free parameter(s).
  ! PAR : REAL(KIND=8)
  !     Equation parameters
  !
  ! Output
  ! ----------
  ! F_out  : REAL(KIND=8), ARRAY
  !    Equation or ODE right-hand-side values

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE
  
  !--------------------------------!
  !     STATE-SPACE DIMENSIONS     !
  !--------------------------------!
  ! State-space dimension of the Yamada model
  INTEGER, PARAMETER          :: xdim = 3
  ! Parameter-space dimension of the Yamada model
  INTEGER, PARAMETER          :: pdim = 4

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! Array indicating the free parameter(s)
  INTEGER, INTENT(IN)         :: ICP(*)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! Array indicating the free parameter(s)
  INTEGER, INTENT(IN)         :: IJAC

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: F_out(NDIM)

  ! State space Jacobian
  REAL(KIND=8), INTENT(INOUT) :: DFDU(NDIM, NDIM)
  ! Parameter space Jacobian
  REAL(KIND=8), INTENT(INOUT) :: DFDP(NDIM, *)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! State space Jacobian
  REAL(KIND=8)                :: J_DFDX(xdim, xdim)
  ! Segment 1: State-space and adjoint-space vectors
  REAL(KIND=8)                :: x1_vec(xdim), w1_vec(xdim)
  ! Segment 2: State-space and adjoint-space vectors
  REAL(KIND=8)                :: x2_vec(xdim), w2_vec(xdim)
  ! Segment 3: State-space vector
  REAL(KIND=8)                :: x3_vec(xdim)
  ! Segment 4: State-space vector
  REAL(KIND=8)                :: x4_vec(xdim)
  ! Yamada model parameters
  REAL(KIND=8)                :: gamma, A_pump, B, a
  ! Phase resetting parameters
  REAL(KIND=8)                :: T, k, theta_old, theta_new

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! Segment 1: State-space and adjoint-space vectors
  x1_vec    = U(1 : xdim)
  w1_vec    = U(xdim+1 : 2*xdim)
  ! Segment 2: State-space and adjoint-space vectors
  x2_vec    = U((2*xdim)+1 : 3*xdim)
  w2_vec    = U((3*xdim)+1 : 4*xdim)
  ! Segment 3: State-space vector
  x3_vec    = U((4*xdim)+1 : 5*xdim)
  ! Segment 4: State-space vector
  x4_vec    = U((5*xdim)+1 : 6*xdim)

  ! Yamada model parameters
  gamma     = PAR(1)
  A_pump    = PAR(2)
  B         = PAR(3)
  a         = PAR(4)

  ! Period of orbit
  T         = PAR(pdim+1)
  ! Integer number of perioids
  k         = PAR(pdim+2)
  ! Phase where the perturbation is applied
  theta_old = PAR(pdim+5)
  ! Phase where segments comes back to \Gamma
  theta_new = PAR(pdim+6)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  !-------------------!
  !     Segment 1     !
  !-------------------!
  ! Vector field
  F_out(1)   = gamma * (A_pump - x1_vec(1) - (x1_vec(1) * x1_vec(3)))
  F_out(2)   = gamma * (B - x1_vec(2) - (a * x1_vec(2) * x1_vec(3)))
  F_out(3)   = (x1_vec(1) - x1_vec(2) - 1.0d0) * x1_vec(3)

  ! Adjoint equation
  F_out(4)   = (-x1_vec(3) * w1_vec(3)) + &
             & (gamma * w1_vec(1) * (x1_vec(3) + 1.0d0))
  F_out(5)   = (x1_vec(3) * w1_vec(3)) + &
             & (gamma * w1_vec(2) * ((a * x1_vec(3)) + 1.0d0))
  F_out(6)   = -w1_vec(3) * (x1_vec(1) - x1_vec(2) - 1.0d0) + &
             & gamma * ((w1_vec(1) * x1_vec(1)) + (w1_vec(2) * a * x1_vec(2)))

  ! Multiply by T \theta_new
  F_out(1:6) = F_out(1:6) * T * theta_new

  !-------------------!
  !     Segment 2     !
  !-------------------!
  ! Vector field
  F_out(7)    = gamma * (A_pump - x2_vec(1) - (x2_vec(1) * x2_vec(3)))
  F_out(8)    = gamma * (B - x2_vec(2) - (a * x2_vec(2) * x2_vec(3)))
  F_out(9)    = (x2_vec(1) - x2_vec(2) - 1.0d0) * x2_vec(3)

  ! Adjoint equation
  F_out(10)   = (-x2_vec(3) * w2_vec(3)) + &
              & (gamma * w2_vec(1) * (x2_vec(3) + 1.0d0))
  F_out(11)   = (x2_vec(3) * w2_vec(3)) + &
              & (gamma * w2_vec(2) * ((a * x2_vec(3)) + 1.0d0))
  F_out(12)   = -w2_vec(3) * (x2_vec(1) - x2_vec(2) - 1.0d0) + &
              & gamma * ((w2_vec(1) * x2_vec(1)) + (w2_vec(2) * a * x2_vec(2)))

  ! Multiply by T (1 - \theta_new)
  F_out(7:12) = F_out(7:12) * T * (1.0d0 - theta_new)

  !-------------------!
  !     Segment 3     !
  !-------------------!
  ! Vector field
  F_out(13)    = gamma * (A_pump - x3_vec(1) - (x3_vec(1) * x3_vec(3)))
  F_out(14)    = gamma * (B - x3_vec(2) - (a * x3_vec(2) * x3_vec(3)))
  F_out(15)    = (x3_vec(1) - x3_vec(2) - 1.0d0) * x3_vec(3)

  ! Multiply by T (1 - \theta_old)
  F_out(13:15) = F_out(13:15) * T * (1.0d0 - theta_old)

  !-------------------!
  !     Segment 4     !
  !-------------------!
  ! Vector field
  F_out(16)    = gamma * (A_pump - x4_vec(1) - (x4_vec(1) * x4_vec(3)))
  F_out(17)    = gamma * (B - x4_vec(2) - (a * x4_vec(2) * x4_vec(3)))
  F_out(18)    = x4_vec(3) * (x4_vec(1) - x4_vec(2) - 1.0d0)

  ! Multiply by k T
  F_out(16:18) = F_out(16:18) * k * T

  IF (IJAC .EQ. 0) RETURN
  !============================================================================!
  !                            JACOBIAN: STATE-SPACE                           !
  !============================================================================!


  IF (IJAC .EQ. 1) RETURN
  !============================================================================!
  !                          JACOBIAN: PARAMETER-SPACE                         !
  !============================================================================!

  ! End of subroutine
END SUBROUTINE FUNC

!==============================================================================!
!                             BOUNDARY CONDITIONS                              !
!==============================================================================!
SUBROUTINE BCND(NDIM, PAR, ICP, NBC, U0, U1, FB, IJAC, DBC)
  ! Boundary conditions of the four segments of the phase resetting problem:
  !                        x1(0) - x2(1) = 0 ,
  !                        x1(1) - x2(0) = 0 ,
  !                        e1 . F(x1(0)) = 0 ,
  !                        w1(0) - w2(1) = 0 ,
  !                   mu_s w1(1) - w2(0) = 0 ,
  !                          |w2(0)| - 1 = 0 ,
  !                        x3(1) - x1(0) = 0 ,
  !                x4(0) - x3(0) - A d_r = 0 ,
  !              (x4(1) - x2(0)) . w2(0) = 0 ,
  !             | x4(1) - x2(0) | - \eta = 0 .
  !
  ! Input
  ! ----------
  ! NDIM : INTEGER
  !     Dimension of the algebraic or ODE system
  ! PAR  : REAL(KIND=8), ARRAY
  !     Array of the equation parameters
  ! ICP  : INTEGER, ARRAY
  !     Array indicating the free parameters
  ! NBC  : INTEGER
  !     Number of boundary conditions
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! IJAC : INTEGER
  !     Something for the Jacobian (not usually used)
  !
  ! Output
  ! ----------
  ! FB   : REAL(KIND=8), ARRAY
  !     The value of the boundary condition functions
  ! DBC  : REAL(KIND=8), ARRAY
  !     Something else for the Jacobian (not usually used)

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE
  
  !--------------------------------!
  !     STATE-SPACE DIMENSIONS     !
  !--------------------------------!
  ! State-space dimension of the Yamada model
  INTEGER, PARAMETER          :: xdim = 3
  ! Parameter-space dimension of the Yamada model
  INTEGER, PARAMETER          :: pdim = 4

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! Equation parameters.
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! Array indicating the free parameter(s)
  INTEGER, INTENT(IN)         :: ICP(*)
  ! Number of boundary conditions
  INTEGER, INTENT(IN)         :: NBC
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(NDIM)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(NDIM)
  ! Something for the Jacobian
  INTEGER, INTENT(IN)         :: IJAC

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the boundary condition functions
  REAL(KIND=8), INTENT(OUT)   :: FB(NBC)
  ! Something else for the Jacobian (not really used)
  REAL(KIND=8), INTENT(INOUT) :: DBC(NBC, *)

  !---------------------------!
  !     SUBROUTINE THINGS     !
  !---------------------------!
  ! Segment 1: Initial and final vectors
  REAL(KIND=8)                :: x0_seg1(xdim), w0_seg1(xdim)
  REAL(KIND=8)                :: x1_seg1(xdim), w1_seg1(xdim)
  ! Segment 2: Initial and final vectors
  REAL(KIND=8)                :: x0_seg2(xdim), w0_seg2(xdim)
  REAL(KIND=8)                :: x1_seg2(xdim), w1_seg2(xdim)
  ! Segment 3: Initial and final vectors
  REAL(KIND=8)                :: x0_seg3(xdim)
  REAL(KIND=8)                :: x1_seg3(xdim)
  ! Segment 4: Initial and final vectors
  REAL(KIND=8)                :: x0_seg4(xdim)
  REAL(KIND=8)                :: x1_seg4(xdim)
  ! Yamada model parameters
  REAL(KIND=8)                :: gamma, A_pump, B, a
  ! Parameters
  REAL(KIND=8)                :: mu_s, eta
  ! Parameters
  REAL(KIND=8)                :: A_perturb, theta_perturb, phi_perturb
  ! Directional vector
  REAL(KIND=8)                :: d_perturb(xdim)

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  !-------------------------!
  !     Initial Vectors     !
  !-------------------------!
  ! Segment 1
  x0_seg1 = U0(1 : xdim)
  w0_seg1 = U0(xdim+1 : 2*xdim)

  ! Segment 2
  x0_seg2 = U0(2*xdim+1 : 3*xdim)
  w0_seg2 = U0(3*xdim+1 : 4*xdim)

  ! Segment 3
  x0_seg3 = U0(4*xdim+1 : 5*xdim)

  ! Segment 4
  x0_seg4 = U0(5*xdim+1 : 6*xdim)

  !-----------------------!
  !     Final Vectors     !
  !-----------------------!
  ! Segment 1
  x1_seg1 = U1(1 : xdim)
  w1_seg1 = U1(xdim+1 : 2*xdim)

  ! Segment 2
  x1_seg2 = U1(2*xdim+1 : 3*xdim)
  w1_seg2 = U1(3*xdim+1 : 4*xdim)

  ! Segment 3
  x1_seg3 = U1(4*xdim+1 : 5*xdim)

  ! Segment 4
  x1_seg4 = U1(5*xdim+1 : 6*xdim)

  !--------------------!
  !     Parameters     !
  !--------------------!
  ! Yamada model parameters
  gamma         = PAR(1)
  A_pump        = PAR(2)
  B             = PAR(3)
  a             = PAR(4)

  ! Stable Floquet eigenvalue
  mu_s          = PAR(pdim+3)
  ! Distance from perturbed segment to \Gamma
  eta           = PAR(pdim+4)
  ! Perturbation amplitude
  A_perturb     = PAR(pdim+7)
  ! Angle at which perturbation is applied
  theta_perturb = PAR(pdim+8)
  ! Azimuthal angle at which perturbation is applied
  phi_perturb   = PAR(pdim+9)

  ! Perturbation directional vector
  d_perturb(1) = COS(theta_perturb)
  d_perturb(2) = 0.0d0
  d_perturb(3) = SIN(theta_perturb)

  !============================================================================!
  !                        BOUNDARY CONDITION ENCODING                         !
  !============================================================================!
  !---------------------------------!
  !     Segment 1 and Segment 2     !
  !---------------------------------!
  ! Boundary conditions - Segments 1 and 2
  FB(1:3)   = x0_seg1 - x1_seg2
  FB(4:6)   = x1_seg1 - x0_seg2
  FB(7)     = gamma * (A_pump - x0_seg1(1) - (x0_seg1(1) * x0_seg1(3)))

  ! Adjoint Boundary Conditions - Segments 1 and 2
  FB(8:10)  = w0_seg1 - w1_seg2
  FB(11:13) = (mu_s * w0_seg2) - w1_seg1
  FB(14)    = NORM2(w0_seg2) - 1.0d0

  !-------------------!
  !     Segment 3     !
  !-------------------!
  ! Boundary conditions - Segments 1 and 3
  FB(15:17) = x1_seg3 - x0_seg1

  !-------------------!
  !     Segment 4     !
  !-------------------!
  FB(18:20) = x0_seg4 - x0_seg3 - (A_perturb * d_perturb)
  FB(21)    = DOT_PRODUCT(x1_seg4 - x0_seg2, w0_seg2)
  FB(22)    = NORM2(x1_seg4 - x0_seg2) - eta

END SUBROUTINE BCND

!==============================================================================!
!                              INITIAL CONDITIONS                              !
!==============================================================================!
SUBROUTINE STPNT(NDIM, U, PAR, T)
  ! Vector field encoding of the Yamada model of a saturable
  ! laser. The vector field is described by the set of
  ! coupled ODEs:
  !          G' = \gamma (A - G - G I) ,
  !          Q' = \gamma (B - Q - a Q I) ,
  !          I' = (G - Q - 1) I .
  !
  ! Input
  ! ----------
  ! NDIM    : INTEGER
  !     Dimension of the algebraic or ODE system
  ! T       : REAL(KIND=8)
  !     Time array for time- or space-dependent fields
  !
  ! Output
  ! ----------
  ! U   : REAL(KIND=8)
  !     State variables of the vector field
  ! PAR : REAL(KIND=8)
  !     Equation parameters

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! Time value for time-dependent solutions
  REAL(KIND=8), INTENT(IN)    :: T

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! State variables of the vector field
  REAL(KIND=8), INTENT(INOUT)    :: U(NDIM)
  ! Equation parameters.
  REAL(KIND=8), INTENT(INOUT)    :: PAR(*)

  !============================================================================!
  !                          INITIAL FIELD COMPONENTS                          !
  !============================================================================!

  ! End of subroutine
END SUBROUTINE STPNT

!==============================================================================!
!                             INTEGRAL CONDITIONS                              !
!==============================================================================!
SUBROUTINE ICND(NDIM, PAR, ICP, NINT, U, UOLD, UDOT, UPOLD, FI, IJAC, DINT)
  ! Integral conditions...
  !
  ! Input
  ! ----------
  ! NDIM  : INTEGER
  !     Dimension of the algebraic or ODE system.
  ! PAR   : REAL(KIND=8), ARRAY
  !     Array of the equation parameters
  ! ICP   : INTEGER, ARRAY
  !     Array indicating the free parameters
  ! NINT  : INTEGER
  !     Number of integral conditions
  ! U     : REAL(KIND=8), ARRAY
  !     Value of the vector function U at 'time' t
  !
  ! The following input arguments, which are normally not needed,
  ! correspond to the preceding point on the solution branch
  ! UOLD  : REAL(KIND=8), ARRAY
  !     The state vector at 'time' t
  ! UDOT  : REAL(KIND=8), ARRAY
  !     Derivate of UOLD w.r.t arclength
  ! UPOLD : REAL(KIND=8), ARRAY
  !     Derivate of UOLD w.r.t 'time'
  ! IJAC  : INTEGER
  !     Something for the Jacobian (not usually used)
  ! 
  ! Output
  ! ----------
  ! FI    : REAL(KIND=8), ARRAY
  !     The value of the vector integrand
  ! DINT  : REAL(KIND=8), ARRAY
  !     Something else for the Jacobian (not usually used)

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! Equation parameters.
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! Array indicating the free parameter(s)
  INTEGER, INTENT(IN)         :: ICP(*)
  ! Number of boundary conditions
  INTEGER, INTENT(IN)         :: NINT
  ! Value of the vector function U at 'time' t
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! The state vector at 'time' t
  REAL(KIND=8), INTENT(IN)    :: UOLD(NDIM)
  ! Derivate of UOLD w.r.t arclength
  REAL(KIND=8), INTENT(IN)    :: UDOT(NDIM)
  ! Derivate of UOLD w.r.t 'time'
  REAL(KIND=8), INTENT(IN)    :: UPOLD(NDIM)
  ! Something for the Jacobian
  INTEGER, INTENT(IN)         :: IJAC

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the vector integrand
  REAL(KIND=8), INTENT(OUT)   :: FI(NINT)
  ! Something else for the Jacobian (not really used)
  REAL(KIND=8), INTENT(INOUT) :: DINT(NINT, *)

  !============================================================================!
  !                        INTEGRAL CONDITION ENCODING                         !
  !============================================================================!
  ! Integral condition encoding
  ! FI(1) = ...
  ! FI(2) = ...

END SUBROUTINE ICND

!==============================================================================!
!                              ALGEBRAIC PROBLEMS                              !
!==============================================================================!
SUBROUTINE FOPT(NDIM, U, ICP, PAR, IJAC, FS, DFDU, DFDP)
  ! Defines the objective function for algebraic optimization problems
  !
  ! Supplied variables :
  !      NDIM   :   Dimension of the state equation
  !      U      :   The state vector
  !      ICP    :   Indices of the control parameters
  !      PAR    :   The vector of control parameters
  !
  ! Values to be returned :
  !      FS      :   The value of the objective function
  !
  ! Normally unused Jacobian argument : IJAC, DFDP

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! The state vector
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! Indices of the control parameters
  INTEGER, INTENT(IN)         :: ICP(*)
  ! The vector of control parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! Something for the Jacobian
  INTEGER, INTENT(IN)         :: IJAC

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the vector integrand
  REAL(KIND=8), INTENT(OUT)   :: FS
  ! Something else for the Jacobian (not really used)
  REAL(KIND=8), INTENT(INOUT) :: DFDU(NDIM)
  REAL(KIND=8), INTENT(INOUT) :: DFDP(*)

  !X FS=

END SUBROUTINE FOPT

SUBROUTINE PVLS(NDIM, U, PAR)
  ! NOTE : 
  ! Parameters set in this subroutine should be considered as ``solution 
  ! measures'' and be used for output purposes only.
  ! 
  ! They should never be used as `true'' continuation parameters. 
  !
  ! They may, however, be added as ``over-specified parameters'' in the 
  ! parameter list associated with the AUTO-Constant NICP, in order to 
  ! print their values on the screen and in the ``p.xxx file.
  !
  ! They may also appear in the list associated with AUTO-Constant NUZR.
  !
  !---------------------------------------------------------------------- 
  ! For algebraic problems the argument U is, as usual, the state vector.
  ! For differential equations the argument U represents the approximate 
  ! solution on the entire interval [0,1]. In this case its values must 
  ! be accessed indirectly by calls to GETP, as illustrated below.
  !---------------------------------------------------------------------- 
  !
  ! Set PAR(2) equal to the L2-norm of U(1)
  !X PAR(2)=GETP('NRM',1,U)
  !
  ! Set PAR(3) equal to the minimum of U(2)
  !X PAR(3)=GETP('MIN',2,U)
  !
  ! Set PAR(4) equal to the value of U(2) at the left boundary.
  !X PAR(4)=GETP('BV0',2,U)
  !
  ! Set PAR(5) equal to the pseudo-arclength step size used.
  !X PAR(5)=GETP('STP',1,U)
  !
  !---------------------------------------------------------------------- 
  ! The first argument of GETP may be one of the following:
  !        'NRM' (L2-norm),     'MAX' (maximum),
  !        'INT' (integral),    'BV0 (left boundary value),
  !        'MIN' (minimum),     'BV1' (right boundary value).
  !
  ! Also available are
  !   'STP' (Pseudo-arclength step size used).
  !   'FLD' (`Fold function', which vanishes at folds).
  !   'BIF' (`Bifurcation function', which vanishes at singular points).
  !   'HBF' (`Hopf function'; which vanishes at Hopf points).
  !   'SPB' ( Function which vanishes at secondary periodic bifurcations).
  !---------------------------------------------------------------------- 
  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! Equation parameters.
  REAL(KIND=8), INTENT(INOUT) :: PAR(*)

  
END SUBROUTINE PVLS