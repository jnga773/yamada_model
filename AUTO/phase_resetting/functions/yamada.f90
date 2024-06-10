!==============================================================================!
!                             AUTO-07: YAMADA MODEL                            !
!==============================================================================!
! This module file contains the subroutines neeed to calculate things for the
! Yamada model using AUTO-07p.
!
! This module contains the following subroutines:
! - FIELD      : Vector field encoding of the Yamada model.
! - FIELD_DFDX : Encoding of the state-space Jacobian.
! - FIELD_DFDP : Encoding of the parameter-space Jacobian.
!
! - FUNC       : The encoded vector field and Jacobians of the continuation
!                problem.
!
! The following subroutines are needed for AUTO, but are not used.
! - STPNT      : Define and set initial conditions of the vector field.
! - BCND       : Define the various boundary conditions (if necessary).
! - ICND       : Define the various integral conditions (if necessary).
!
! - FOPT       : Defines the objective function for algebraic optimisation
!                problems (if necessary).
! - PVLS       : Something for algebraic problems (which I don't need to worry
!                about).
!
! The U-vector components are as follows:
!         U(1:3) - State-space components of the vector field (equilibrium
!                  point or periodic orbit, depending on the calculation).
!
! The parameters (PAR) are as follows:
!       PAR(1:4) - Parameters of the Yamada model (gamma, A, B, a),
!       PAR(11)  - (FOR PERIODIC ORBIT CALCULATION ONLY) Period of the segment.

!==============================================================================!
!                     CHANGE THESE VARIABLES AND FUNCTIONS                     !
!==============================================================================!
MODULE VECTOR_FIELD
  ! This module contains the global variables for 'xdim' and 'pdim' (in COCO
  ! language), the original dimensions of the state vector and parameter vector
  ! for the default vector field.

  IMPLICIT NONE

  ! Original dimension of vector field
  INTEGER, PARAMETER :: xdim = 3
  ! Original dimension of parameter vector
  INTEGER, PARAMETER :: pdim = 4

  CONTAINS

  SUBROUTINE FIELD(x_in, p_in, F_out)
    ! Vector field encoding of the Yamada model of a saturable
    ! laser. The vector field is described by the set of
    ! coupled ODEs:
    !          G' = \gamma (A - G - G I) ,
    !          Q' = \gamma (B - Q - a Q I) ,
    !          I' = (G - Q - 1) I .
    !
    ! Input
    ! ----------
    ! x_in   : REAL(KIND=8), ARRAY
    !     State variables of the vector field
    ! par_in : REAL(KIND=8), ARRAY
    !     Equation parameters.
    !
    ! Output
    ! ----------
    ! F_out  : REAL(KIND=8), ARRAY
    !    Equation or ODE right-hand-side values

    !============================================================================!
    !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
    !============================================================================!

    IMPLICIT NONE

    !-------------------------!
    !     INPUT ARGUMENTS     !
    !-------------------------!
    ! State variables of the vector field
    REAL(KIND=8), INTENT(IN)    :: x_in(xdim)
    ! Equation parameters
    REAL(KIND=8), INTENT(IN)    :: p_in(*)

    !--------------------------!
    !     OUTPUT ARGUMENTS     !
    !--------------------------!
    ! Equation or ODE right-hand-side values
    REAL(KIND=8), INTENT(OUT)   :: F_out(xdim)

    !--------------------------!
    !     SUBROUTINE STUFF     !
    !--------------------------!
    ! State space variables
    REAL(KIND=8)                :: G, Q, I
    ! Vector field parameters
    REAL(KIND=8)                :: gamma, A_pump, B, a

    !============================================================================!
    !                            VECTOR FIELD ENCODING                           !
    !============================================================================!
    ! Grab the state-space variables from x_in
    G = x_in(1)
    Q = x_in(2)
    I = x_in(3)

    ! Grab the parameter variables from p_in
    gamma  = p_in(1)
    A_pump = p_in(2)
    B      = p_in(3)
    a      = p_in(4)

    ! The system of equations
    F_out(1) = gamma * (A_pump - G - (G * I))
    F_out(2) = gamma * (B - Q - (a * Q * I))
    F_out(3) = (G - Q - 1.0d0) * I

  END SUBROUTINE FIELD

  SUBROUTINE FIELD_DFDX(x_in, p_in, J_out)
    ! Encoding of the state-space Jacobian of the Yamada model.
    !
    ! Input
    ! ----------
    ! x_in   : REAL(KIND=8), ARRAY
    !     State variables of the vector field
    ! par_in : REAL(KIND=8), ARRAY
    !     Equation parameters.
    !
    ! Output
    ! ----------
    ! J_out  : REAL(KIND=8), ARRAY(NDIM, NDIM)
    !    State-space Jacobian

    !============================================================================!
    !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
    !============================================================================!

    IMPLICIT NONE

    !-------------------------!
    !     INPUT ARGUMENTS     !
    !-------------------------!
    ! State variables of the vector field
    REAL(KIND=8), INTENT(IN)    :: x_in(xdim)
    ! Equation parameters
    REAL(KIND=8), INTENT(IN)    :: p_in(*)

    !--------------------------!
    !     OUTPUT ARGUMENTS     !
    !--------------------------!
    ! Equation or ODE right-hand-side values
    REAL(KIND=8), INTENT(OUT)   :: J_out(xdim, xdim)

    !--------------------------!
    !     SUBROUTINE STUFF     !
    !--------------------------!
    ! State space variables
    REAL(KIND=8)                :: G, Q, I
    ! Vector field parameters
    REAL(KIND=8)                :: gamma, A_pump, B, a

    !============================================================================!
    !                            VECTOR FIELD ENCODING                           !
    !============================================================================!
    ! Grab the state-space variables from U
    G = x_in(1)
    Q = x_in(2)
    I = x_in(3)

    ! Grab the parameter variables from p_in
    gamma  = p_in(1)
    A_pump = p_in(2)
    B      = p_in(3)
    a      = p_in(4)

    ! State-space Jacobian
    J_out(1, 1) = -gamma * (1.0d0 + I)
    J_out(1, 2) = 0.0d0
    J_out(1, 3) = -gamma * G

    J_out(2, 1) = 0.0d0
    J_out(2, 2) = -gamma * (1.0d0 + (a * I))
    J_out(2, 3) = -gamma * a * Q

    J_out(3, 1) = I
    J_out(3, 2) = -I
    J_out(3, 3) = G - Q - 1.0d0

  END SUBROUTINE FIELD_DFDX

  SUBROUTINE FIELD_DFDP(x_in, p_in, J_out)
    ! Encoding of the parameter-space Jacobian of the Yamada model.
    !
    ! Input
    ! ----------
    ! x_in   : REAL(KIND=8), ARRAY
    !     State variables of the vector field
    ! par_in : REAL(KIND=8), ARRAY
    !     Equation parameters.
    !
    ! Output
    ! ----------
    ! J_out  : REAL(KIND=8), ARRAY(NDIM, NDIM)
    !    Parameter-space Jacobian

    !============================================================================!
    !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
    !============================================================================!

    IMPLICIT NONE

    !-------------------------!
    !     INPUT ARGUMENTS     !
    !-------------------------!
    ! State variables of the vector field
    REAL(KIND=8), INTENT(IN)    :: x_in(xdim)
    ! Equation parameters
    REAL(KIND=8), INTENT(IN)    :: p_in(*)

    !--------------------------!
    !     OUTPUT ARGUMENTS     !
    !--------------------------!
    ! Equation or ODE right-hand-side values
    REAL(KIND=8), INTENT(OUT)   :: J_out(xdim, *)

    !--------------------------!
    !     SUBROUTINE STUFF     !
    !--------------------------!
    ! State space variables
    REAL(KIND=8)                :: G, Q, I
    ! Vector field parameters
    REAL(KIND=8)                :: gamma, A_pump, B, a

    !============================================================================!
    !                            VECTOR FIELD ENCODING                           !
    !============================================================================!
    ! Grab the state-space variables from U
    G = x_in(1)
    Q = x_in(2)
    I = x_in(3)

    ! Grab the parameter variables from p_in
    gamma  = p_in(1)
    A_pump = p_in(2)
    B      = p_in(3)
    a      = p_in(4)

    ! Parameter-space Jacobian
    J_out(1, 1) = A_pump - G - (G * I)
    J_out(1, 2) = gamma
    J_out(1, 3) = 0.0d0
    J_out(1, 4) = 0.0d0

    J_out(2, 1) = B - Q - (a * Q * I)
    J_out(2, 2) = 0.0d0
    J_out(2, 3) = gamma
    J_out(2, 4) = -gamma * Q * I

    J_out(3, 1) = 0.0d0
    J_out(3, 2) = 0.0d0
    J_out(3, 3) = 0.0d0
    J_out(3, 4) = 0.0d0

  END SUBROUTINE FIELD_DFDP
END MODULE VECTOR_FIELD

!==============================================================================!
!                          VECTOR FIELD AND JACOBIANS                          !
!==============================================================================!
SUBROUTINE FUNC(NDIM, U, ICP, PAR, IJAC, F_out, DFDU, DFDP)
  ! Vector field encoding of continuation problem.
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
  ! Import global variables
  USE VECTOR_FIELD

  IMPLICIT NONE

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

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  ! Get Yamada Encoding
  CALL FIELD(U(1:xdim), PAR, F_out(1:xdim))

  ! State space Jacobian
  CALL FIELD_DFDX(U(1:xdim), PAR, J_DFDX(:, :))

  !============================================================================!
  !                            HARD-CODED JACOBIANS                            !
  !============================================================================!
  IF (IJAC .EQ. 0) RETURN

  ! If NDIM = 3, calculate base Jacobians
  IF (NDIM .EQ. xdim) THEN

    ! State-space Jacobian
    CALL FIELD_DFDX(U, PAR, DFDU)

    IF (IJAC .EQ. 1) RETURN

    ! Parameter-space Jacobian
    CALL FIELD_DFDP(U, PAR, DFDP)
  END IF

  ! End of subroutine
END SUBROUTINE FUNC

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
!                              BOUNDARY CONDITIONS                             !
!==============================================================================!
SUBROUTINE BCND(NDIM, PAR, ICP, NBC, U0, U1, FB, IJAC, DBC)
  ! Boundary conditions...
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

  !============================================================================!
  !                        BOUNDARY CONDITION ENCODING                         !
  !============================================================================!
  ! Boundary condition encoding
  ! FB(1) = ...
  ! FB(2) = ...

END SUBROUTINE BCND

!==============================================================================!
!                              INITIAL CONDITIONS                              !
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