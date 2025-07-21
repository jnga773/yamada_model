!==============================================================================!
!                           AUTO-07: PHASE RESETTING                           !
!==============================================================================!
! This module file contains the subroutines neeed to calculate the phase
! resetting problem of the Yamada model in AUTO.
!
! This module contains the following subroutines:
! - FIELD          : Vector field encoding of the Yamada model.
! - FIELD_DFDX     : Encoding of the state-space Jacobian.
! - FIELD_DFDP     : Encoding of the parameter-space Jacobian.
! - ADJOINT_DFDX   : Encoding of the state-space Jacobian of the adjoint
!                    equations.
!
! - FUNC           : The encoded vector field of the continuation problem.
! - PR_SEGS        : The encoded vector field of the four segments of the
!                    phase resetting problem.
! - PR_SEGS_DFDX   : The encoded state-space Jacobian of the four segments
!                    of the phase resetting problem.
! - PR_SEGS_DFDP   : The encoded parameter-space Jacobian of the four segments
!                    of the phase resetting problem.
!
! - BCND           : Define the various boundary conditions (if necessary).
! - BCS_SEGS       : Encoding of the boundary conditions of the four segments
!                    of the phase resetting problem.
! - BCS_SEGS_DFDU  : Encoding of the state-space Jacobian of the boundary
!                    conditions of the four segments of the phase resetting
!                    problem.
! - BCS_SEGS_DFDP  : Encoding of the parameter-space Jacobian of the boundary
!                    conditions of the four segments of the phase resetting
!                    problem.
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
    ! PAR : REAL(KIND=8), ARRAY
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
    !                              INPUT PARAMETERS                              !
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

    !============================================================================!
    !                            VECTOR FIELD ENCODING                           !
    !============================================================================!
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
    ! PAR : REAL(KIND=8), ARRAY
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
    !                              INPUT PARAMETERS                              !
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

    !============================================================================!
    !                             JACOBIAN ENCODING                              !
    !============================================================================!
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
    ! PAR : REAL(KIND=8), ARRAY
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
    !                              INPUT PARAMETERS                              !
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

    !============================================================================!
    !                             JACOBIAN ENCODING                              !
    !============================================================================!
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

  SUBROUTINE ADJOINT_DFDX(U_in, p_in, J_out)
    ! Encoding of the state-space Jacobian of the adjoint equations.
    !
    ! Input
    ! ----------
    ! U_in   : REAL(KIND=8), ARRAY
    !     State variables of the vector field
    ! PAR : REAL(KIND=8), ARRAY
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
    REAL(KIND=8), INTENT(IN)    :: U_in(2*xdim)
    ! Equation parameters
    REAL(KIND=8), INTENT(IN)    :: p_in(*)

    !--------------------------!
    !     OUTPUT ARGUMENTS     !
    !--------------------------!
    ! Equation or ODE right-hand-side values
    REAL(KIND=8), INTENT(OUT)   :: J_out(xdim, 2*xdim)

    !--------------------------!
    !     SUBROUTINE STUFF     !
    !--------------------------!
    ! State space variables
    REAL(KIND=8)                :: G, Q, I
    ! Adjoint space variables
    REAL(KIND=8)                :: w1, w2, w3
    ! Vector field parameters
    REAL(KIND=8)                :: gamma, A_pump, B, a

    !============================================================================!
    !                              INPUT PARAMETERS                              !
    !============================================================================!
    ! Grab the state-space variables from U_in
    G  = U_in(1)
    Q  = U_in(2)
    I  = U_in(3)
    ! Grab the adjoint-space variables from v
    w1 = U_in(4)
    w2 = U_in(5)
    w3 = U_in(6)

    ! Grab the parameter variables from p_in
    gamma  = p_in(1)
    A_pump = p_in(2)
    B      = p_in(3)
    a      = p_in(4)

    !============================================================================!
    !                             JACOBIAN ENCODING                              !
    !============================================================================!
    ! State-space Jacobian
    J_out(1, 3) = (gamma * w1) + (gamma * a * w2)
    J_out(1, 4) = gamma * (1.0d0 + I)
    J_out(1, 6) = -I
    
    J_out(2, 3) = (gamma * a * w2) + w3
    J_out(2, 5) = gamma * ((a * I) + 1.0d0)
    J_out(2, 6) = I

    J_out(3, 1) = (gamma * a * w2) - w3
    J_out(3, 2) = (gamma * a * w2) + w3
    J_out(3, 4) = gamma * G
    J_out(3, 5) = gamma * a * Q
    J_out(3, 6) = -(G - Q - 1.0d0)
  
  END SUBROUTINE ADJOINT_DFDX

END MODULE VECTOR_FIELD

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

  !============================================================================!
  !                       CONTINUATION PROBLEM ENCODING                        !
  !============================================================================!
  !----------------------!
  !     VECTOR FIELD     !
  !----------------------!
  CALL PR_SEGS(NDIM, U, PAR, F_out)

  IF (IJAC .EQ. 0) RETURN
  !-------------------------------!
  !     JACOBIAN: STATE-SPACE     !
  !-------------------------------!
  CALL PR_SEGS_DFDX(NDIM, U, PAR, DFDU)

  IF (IJAC .EQ. 1) RETURN
  !-----------------------------------!
  !     JACOBIAN: PARAMETER-SPACE     !
  !-----------------------------------!
  CALL PR_SEGS_DFDP(NDIM, U, PAR, DFDP)

  ! End of subroutine
END SUBROUTINE FUNC

SUBROUTINE PR_SEGS(NDIM, U, PAR, F_out)
  ! Vector field encoding for the four segments of the phase resetting problem.
  !
  ! Input
  ! ----------
  ! NDIM   : INTEGER
  !     Dimension of the algebraic or ODE system.
  ! U   : REAL(KIND=8)
  !     State variables of the vector field.
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
  ! Import global variables and vector fields
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: F_out(NDIM)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! State space Jacobian
  REAL(KIND=8)                :: J_DFDX(xdim, xdim)
  ! Segment 1: State-space and adjoint-space vectors
  REAL(KIND=8)                :: x_seg1(xdim), w_seg1(xdim)
  ! Segment 2: State-space and adjoint-space vectors
  REAL(KIND=8)                :: x_seg2(xdim), w_seg2(xdim)
  ! Segment 3: State-space vector
  REAL(KIND=8)                :: x_seg3(xdim)
  ! Segment 4: State-space vector
  REAL(KIND=8)                :: x_seg4(xdim)
  ! Yamada model parameters
  REAL(KIND=8)                :: gamma, A_pump, B, a
  ! Phase resetting parameters
  REAL(KIND=8)                :: T, k, theta_old, theta_new

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! Segment 1: State-space and adjoint-space vectors
  x_seg1    = U(1 : xdim)
  w_seg1    = U(xdim+1 : 2*xdim)
  ! Segment 2: State-space and adjoint-space vectors
  x_seg2    = U((2*xdim)+1 : 3*xdim)
  w_seg2    = U((3*xdim)+1 : 4*xdim)
  ! Segment 3: State-space vector
  x_seg3    = U((4*xdim)+1 : 5*xdim)
  ! Segment 4: State-space vector
  x_seg4    = U((5*xdim)+1 : 6*xdim)

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
  theta_old = PAR(pdim+3)
  ! Phase where segments comes back to \Gamma
  theta_new = PAR(pdim+4)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  !-------------------!
  !     Segment 1     !
  !-------------------!
  ! Calculate vector field
  CALL FIELD(x_seg1, PAR, vec_field)
  ! Calculate Jacobian
  CALL FIELD_DFDX(x_seg1, PAR, J_DFDX)

  ! Vector field
  F_out(1 : xdim) = T * theta_new * vec_field
  ! Adjoint equation
  F_out(xdim+1 : 2*xdim) = -T * theta_new * MATMUL(TRANSPOSE(J_DFDX), w_seg1)

  !-------------------!
  !     Segment 2     !
  !-------------------!
  ! Calculate vector field
  CALL FIELD(x_seg2, PAR, vec_field)
  ! Calculate Jacobian
  CALL FIELD_DFDX(x_seg2, PAR, J_DFDX)

  ! Vector field
  F_out(2*xdim+1 : 3*xdim) = T * (1.0d0 - theta_new) * vec_field
  ! Adjoint equation
  F_out(3*xdim+1 : 4*xdim) = -T * (1.0d0 - theta_new) * MATMUL(TRANSPOSE(J_DFDX), w_seg2)

  !-------------------!
  !     Segment 3     !
  !-------------------!
  ! Calculate vector field
  CALL FIELD(x_seg3, PAR, vec_field)

  ! Vector field
  F_out(4*xdim+1 : 5*xdim) = T * (1.0d0 - theta_old) * vec_field

  !-------------------!
  !     Segment 4     !
  !-------------------!
  ! Calculate vector field
  CALL FIELD(x_seg4, PAR, vec_field)

  ! Vector field
  F_out(5*xdim+1 : 6*xdim) = k * T * vec_field

END SUBROUTINE PR_SEGS

SUBROUTINE PR_SEGS_DFDX(NDIM, U, PAR, J_out)
  ! Encoding of the state-sapce Jacobian for the four segments of the
  ! phase resetting problem.
  !
  ! Input
  ! ----------
  ! NDIM   : INTEGER
  !     Dimension of the algebraic or ODE system.
  ! U   : REAL(KIND=8)
  !     State variables of the vector field.
  ! PAR : REAL(KIND=8)
  !     Equation parameters
  !
  ! Output
  ! ----------
  ! J_out  : REAL(KIND=8), ARRAY
  !    Equation or ODE right-hand-side values

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!
  ! Import global variables and vector fields
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: J_out(NDIM, NDIM)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! State space Jacobian
  REAL(KIND=8)                :: Jacobian(xdim, xdim)
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
  ! Cheaty MATLAB constant things
  REAL(KIND=8)                :: t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12
  REAL(KIND=8)                :: t13, t14, t15, t16, t17, t18, t19, t20, t21
  REAL(KIND=8)                :: t22, t23, t24, t25, t26, t27, t28, t29, t30, t31
  REAL(KIND=8)                :: t32, t33, t34, t35, t36, t37, t38, t39, t40
  REAL(KIND=8)                :: t41, t42, t43

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
  theta_old = PAR(pdim+3)
  ! Phase where segments comes back to \Gamma
  theta_new = PAR(pdim+4)

  !============================================================================!
  !                        CHEATY MATLAB CONSTANT THINGS                       !
  !============================================================================!
  t2  = a * x1_vec(3)
  t3  = a * x2_vec(3)
  t4  = x1_vec(3) + 1.0d0
  t5  = x2_vec(3) + 1.0d0
  t6  = T * k * x4_vec(3)
  t7  = T * theta_new * w1_vec(3)
  t8  = T * theta_new * x1_vec(3)
  t9  = theta_old - 1.0d0
  t10 = theta_new - 1.0d0
  t11 = -x1_vec(1)
  t12 = -x2_vec(1)
  t13 = T * gamma * theta_new * w1_vec(1)
  t14 = T * gamma * theta_new * x1_vec(1)
  t15 = t2 + 1.0d0
  t16 = t3 + 1.0d0
  t17 = T * a * gamma * theta_new * w1_vec(2)
  t18 = T * a * gamma * theta_new * x1_vec(2)
  t19 = -t7
  t20 = -t8
  t21 = T * gamma * t4 * theta_new
  t22 = t11 + x1_vec(2) + 1.0d0
  t23 = t12 + x2_vec(2) + 1.0d0
  t24 = T * t10 * w2_vec(3)
  t25 = T * t9 * x3_vec(3)
  t26 = T * t10 * x2_vec(3)
  t27 = T * gamma * t10 * w2_vec(1)
  t28 = T * gamma * t10 * x2_vec(1)
  t29 = T * gamma * t15 * theta_new
  t30 = T * a * gamma * t10 * w2_vec(2)
  t31 = T * a * gamma * t10 * x2_vec(2)
  t32 = -t24
  t33 = -t26
  t34 = -t27
  t35 = T * gamma * t5 * t10
  t36 = T * t22 * theta_new
  t37 = t7 + t17
  t38 = -t30
  t39 = t13 + t19
  t40 = T * gamma * t10 * t16
  t41 = T * t10 * t23
  t42 = t24 + t34
  t43 = t32 + t38

  !============================================================================!
  !                             JACOBIAN ENCODING                              !
  !============================================================================!
  !-------------------!
  !     Segment 1     !
  !-------------------!
  ! Call state-space Jacobian
  CALL FIELD_DFDX(x1_vec, PAR, Jacobian)
  J_out(1:xdim, 1:xdim) = T * theta_new * Jacobian

  ! Call adjoint-space Jacobian
  J_out(4, 3)   = t39
  J_out(4, 4)   = t21
  J_out(4, 6)   = t20
  J_out(5, 3)   = t37
  J_out(5, 5)   = t29
  J_out(5, 6)   = t8
  J_out(6, 1)   = t39
  J_out(6, 2)   = t37
  J_out(6, 4)   = t14
  J_out(6, 5)   = t18
  J_out(6, 6)   = t36

  !-------------------!
  !     Segment 2     !
  !-------------------!
  ! Call state-space Jacobian
  CALL FIELD_DFDX(x2_vec, PAR, Jacobian)
  J_out(2*xdim+1:3*xdim, 2*xdim+1:3*xdim) = T * (1.0d0 - theta_new) * Jacobian

  ! Call adjoint-space Jacobian
  J_out(10, 9)  = t42
  J_out(10, 10) = -t35
  J_out(10, 12) = t26
  J_out(11, 9)  = t43
  J_out(11, 11) = -t40
  J_out(11, 12) = t33
  J_out(12, 7)  = t42
  J_out(12, 8)  = t43
  J_out(12, 10) = T * gamma * t10 * t12
  J_out(12, 11) = -t31
  J_out(12, 12) = -t41

  !-------------------!
  !     Segment 3     !
  !-------------------!
  ! Call state-space Jacobian
  CALL FIELD_DFDX(x3_vec, PAR, Jacobian)
  J_out(4*xdim+1:5*xdim, 4*xdim+1:5*xdim) = T * (1.0d0 - theta_old) * Jacobian

  !-------------------!
  !     Segment 4     !
  !-------------------!
  ! Call state-space Jacobian
  CALL FIELD_DFDX(x4_vec, PAR, Jacobian)
  J_out(5*xdim+1:6*xdim, 5*xdim+1:6*xdim) = k * T * Jacobian

END SUBROUTINE PR_SEGS_DFDX

SUBROUTINE PR_SEGS_DFDP(NDIM, U, PAR, J_out)
  ! Encoding of the parameter-sapce Jacobian for the four segments of the
  ! phase resetting problem.
  !
  ! Input
  ! ----------
  ! NDIM   : INTEGER
  !     Dimension of the algebraic or ODE system.
  ! U   : REAL(KIND=8)
  !     State variables of the vector field.
  ! PAR : REAL(KIND=8)
  !     Equation parameters
  !
  ! Output
  ! ----------
  ! J_out  : REAL(KIND=8), ARRAY
  !    Equation or ODE right-hand-side values

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!
  ! Import global variables and vector fields
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of the algebraic or ODE system
  INTEGER, INTENT(IN)         :: NDIM
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(NDIM)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: J_out(NDIM, *)

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
  ! Cheaty MATLAB constant things
  REAL(KIND=8)                :: t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12
  REAL(KIND=8)                :: t13, t14, t15, t16, t17, t18, t19, t20, t21
  REAL(KIND=8)                :: t22, t23, t24, t25, t26, t27, t28, t29, t30, t31
  REAL(KIND=8)                :: t32, t33, t34, t35, t36, t37, t38, t39, t40
  REAL(KIND=8)                :: t41, t42, t43, t44, t45

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
  theta_old = PAR(pdim+3)
  ! Phase where segments comes back to \Gamma
  theta_new = PAR(pdim+4)
  
  !============================================================================!
  !                        CHEATY MATLAB CONSTANT THINGS                       !
  !============================================================================!
  t2 = a * x1_vec(3)
  t3 = a * x2_vec(3)
  t4 = x1_vec(1) * x1_vec(3)
  t5 = x2_vec(1) * x2_vec(3)
  t6 = x3_vec(1) * x3_vec(3)
  t7 = x4_vec(1) * x4_vec(3)
  t8 = x1_vec(3) + 1.0d0
  t9 = x2_vec(3) + 1.0d0
  t10 = T * gamma * k
  t11 = T * gamma * theta_new
  t12 = -A_pump
  t13 = -B
  t14 = T * w1_vec(3) * x1_vec(3)
  t15 = T * w2_vec(3) * x2_vec(3)
  t16 = t2 * x1_vec(2)
  t17 = t3 * x2_vec(2)
  t18 = a * x3_vec(2) * x3_vec(3)
  t19 = a * x4_vec(2) * x4_vec(3)
  t20 = theta_new * w1_vec(3) * x1_vec(3)
  t21 = theta_old - 1.0d0
  t22 = theta_new - 1.0d0
  t23 = -x1_vec(1)
  t24 = -x2_vec(1)
  t25 = -x3_vec(1)
  t26 = -x4_vec(1)
  t27 = t2 + 1.0d0
  t28 = t3 + 1.0d0
  t29 = t22 * w2_vec(3) * x2_vec(3)
  t30 = t23 + x1_vec(2) + 1.0d0
  t31 = t24 + x2_vec(2) + 1.0d0
  t32 = t25 + x3_vec(2) + 1.0d0
  t33 = t26 + x4_vec(2) + 1.0d0
  t34 = T * gamma * t21
  t35 = T * gamma * t22
  t36 = t4 + t12 + x1_vec(1)
  t37 = t5 + t12 + x2_vec(1)
  t38 = t6 + t12 + x3_vec(1)
  t39 = t7 + t12 + x4_vec(1)
  t40 = -t34
  t41 = -t35
  t42 = t13 + t16 + x1_vec(2)
  t43 = t13 + t17 + x2_vec(2)
  t44 = t13 + t18 + x3_vec(2)
  t45 = t13 + t19 + x4_vec(2)

  !============================================================================!
  !                             JACOBIAN ENCODING                              !
  !============================================================================!
  !-------------------!
  !     Segment 1     !
  !-------------------!
  J_out(1, 1)  = -T * t36 * theta_new
  J_out(1, 2)  = t11
  J_out(1, 5)  = -gamma * t36 * theta_new
  J_out(1, 8)  = -T * gamma * t36
  J_out(2, 1)  = -T * t42 * theta_new
  J_out(2, 3)  = t11
  J_out(2, 4)  = -t11 * x1_vec(2) * x1_vec(3)
  J_out(2, 5)  = -gamma * t42 * theta_new
  J_out(2, 8)  = -T * gamma * t42
  J_out(3, 5)  = -t30 * theta_new * x1_vec(3)
  J_out(3, 8)  = -T * t30 * x1_vec(3)
  J_out(4, 1)  = T * t8 * theta_new * w1_vec(1)
  J_out(4, 5)  = -t20 + &
               & gamma * t8 * theta_new * w1_vec(1)
  J_out(4, 8)  = -t14 + &
               & T * gamma * t8 * w1_vec(1)
  J_out(5, 1)  = T * t27 * theta_new * w1_vec(2)
  J_out(5, 4)  = t11 * w1_vec(2) * x1_vec(3)
  J_out(5, 5)  = t20 + &
               & gamma * t27 * theta_new * w1_vec(2)
  J_out(5, 8)  = t14 + &
               & T * gamma * t27 * w1_vec(2)
  J_out(6, 1)  = T * theta_new * w1_vec(1) * x1_vec(1) + &
               & T * a * theta_new * w1_vec(2) * x1_vec(2)
  J_out(6, 4)  = t11 * w1_vec(2) * x1_vec(2)
  J_out(6, 5)  = t30 * theta_new * w1_vec(3) + &
               & gamma * theta_new * w1_vec(1) * x1_vec(1) + &
               & a * gamma * theta_new * w1_vec(2) * x1_vec(2)
  J_out(6, 8)  = T * t30 * w1_vec(3) + &
               & T * gamma * w1_vec(1) * x1_vec(1) + &
               & T * a * gamma * w1_vec(2) * x1_vec(2)

  !-------------------!
  !     Segment 2     !
  !-------------------!
  J_out(7, 1)  = T * t22 * t37
  J_out(7, 2)  = t41
  J_out(7, 5)  = gamma * t22 * t37
  J_out(7, 8)  = T * gamma * t37
  J_out(8, 1)  = T * t22 * t43
  J_out(8, 3)  = t41
  J_out(8, 4)  = t35 * x2_vec(2) * x2_vec(3)
  J_out(8, 5)  = gamma * t22 * t43
  J_out(8, 8)  = T * gamma * t43
  J_out(9, 5)  = t22 * t31 * x2_vec(3)
  J_out(9, 8)  = T * t31 * x2_vec(3)
  J_out(10, 1) = -T * t9 * t22 * w2_vec(1)
  J_out(10, 5) = t29 - &
               & gamma * t9 * t22 * w2_vec(1)
  J_out(10, 8) = t15 - &
               & T * gamma * t9 * w2_vec(1)
  J_out(11, 1) = -T * t22 * t28 * w2_vec(2)
  J_out(11, 4) = t41 * w2_vec(2) * x2_vec(3)
  J_out(11, 5) = -t29 - &
               & gamma * t22 * t28 * w2_vec(2)
  J_out(11, 8) = -t15 - &
               & T * gamma * t28 * w2_vec(2)
  J_out(12, 1) = T * t22 * t24 * w2_vec(1) - &
               & T * a * t22 * w2_vec(2) * x2_vec(2)
  J_out(12, 4) = t41 * w2_vec(2) * x2_vec(2)
  J_out(12, 5) = -t22 * t31 * w2_vec(3) + &
               & gamma * t22 * t24 * w2_vec(1) - &
               & a * gamma * t22 * w2_vec(2) * x2_vec(2)
  J_out(12, 8) = -T * t31 * w2_vec(3) + &
               & T * gamma * t24 * w2_vec(1) - &
               & T * a * gamma * w2_vec(2) * x2_vec(2)
  
  !-------------------!
  !     Segment 3     !
  !-------------------!
  J_out(13, 1) = T * t21 * t38
  J_out(13, 2) = t40
  J_out(13, 5) = gamma * t21 * t38
  J_out(13, 7) = T * gamma * t38
  J_out(14, 1) = T * t21 * t44
  J_out(14, 3) = t40
  J_out(14, 4) = t34 * x3_vec(2) * x3_vec(3)
  J_out(14, 5) = gamma * t21 * t44
  J_out(14, 7) = T * gamma * t44
  J_out(15, 5) = t21 * t32 * x3_vec(3)
  J_out(15, 7) = T * t32 * x3_vec(3)

  !-------------------!
  !     Segment 4     !
  !-------------------!
  J_out(16, 1) = -T * k * t39
  J_out(16, 2) = t10
  J_out(16, 5) = -gamma * k * t39
  J_out(16, 6) = -T * gamma * t39
  J_out(17, 1) = -T * k * t45
  J_out(17, 3) = t10
  J_out(17, 4) = -t10 * x4_vec(2) * x4_vec(3)
  J_out(17, 5) = -gamma * k * t45
  J_out(17, 6) = -T * gamma * t45
  J_out(18, 5) = -k * t33 * x4_vec(3)
  J_out(18, 6) = -T * t33 * x4_vec(3)

END SUBROUTINE PR_SEGS_DFDP

!==============================================================================!
!                             BOUNDARY CONDITIONS                              !
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
  ! Import global variables and vector fields
  USE VECTOR_FIELD

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
  !-----------------------------!
  !     BOUNDARY CONDITIONS     !
  !-----------------------------!
  CALL BCS_SEGS(NDIM, NBC, U0, U1, PAR, FB)

  IF (IJAC .EQ. 0) RETURN
  !-------------------------------!
  !     JACOBIAN: STATE-SPACE     !
  !-------------------------------!
  ! Jacobian of the boundary condition with respect to the state-space variables
  CALL BCS_SEGS_DFDU(NDIM, NBC, U0, U1, PAR, DBC(1:NBC, 1:NDIM))

  IF (IJAC .EQ. 1) RETURN
  !-----------------------------------!
  !     JACOBIAN: PARAMETER-SPACE     !
  !-----------------------------------!
  ! Jacobian of the boundary condition with respect to the parameters
  CALL BCS_SEGS_DFDP(NDIM, NBC, U0, U1, PAR, DBC(1:NBC, 2*NDIM+1:2*pdim+9))

END SUBROUTINE BCND

SUBROUTINE BCS_SEGS(NDIM, NBC, U0, U1, PAR, BCS_out)
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
  ! NBC  : INTEGER
  !     Number of boundary conditions
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! PAR  : REAL(KIND=8), ARRAY
  !     Array of the equation parameters
  !
  ! Output
  ! ----------
  ! BCS_out : REAL(KIND=8), ARRAY
  !     The value of the boundary condition functions

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
  ! Number of boundary conditions
  INTEGER, INTENT(IN)         :: NBC
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(NDIM)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(NDIM)
  ! Equation parameters.
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the boundary condition functions
  REAL(KIND=8), INTENT(OUT)   :: BCS_out(NBC)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
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
  ! pi
  REAL(KIND=8), PARAMETER    :: pi = 3.14159265358979323846d0
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
  mu_s          = PAR(pdim+5)
  ! Distance from perturbed segment to \Gamma
  eta           = PAR(pdim+6)
  ! Perturbation amplitude
  A_perturb     = PAR(pdim+7)
  ! Angle at which perturbation is applied
  theta_perturb = PAR(pdim+8)
  ! Azimuthal angle at which perturbation is applied
  phi_perturb   = PAR(pdim+9)

  ! Perturbation directional vector
  d_perturb(1) = COS(theta_perturb * (2.0d0 * pi))
  d_perturb(2) = 0.0d0
  d_perturb(3) = SIN(theta_perturb * (2.0d0 * pi))

  !============================================================================!
  !                        BOUNDARY CONDITION ENCODING                         !
  !============================================================================!
  !---------------------------------!
  !     Segment 1 and Segment 2     !
  !---------------------------------!
  ! Boundary conditions - Segments 1 and 2
  BCS_out(1:3)   = x0_seg1 - x1_seg2
  BCS_out(4:6)   = x1_seg1 - x0_seg2
  BCS_out(7)     = gamma * (A_pump - x0_seg1(1) - (x0_seg1(1) * x0_seg1(3)))

  ! Adjoint Boundary Conditions - Segments 1 and 2
  ! BCS_out(8:10)  = w0_seg1 - w1_seg2
  ! BCS_out(11:13) = (mu_s * w0_seg2) - w1_seg1
  ! BCS_out(14)    = NORM2(w0_seg2) - 1.0d0
  BCS_out(8:10)  = w1_seg1 - w0_seg2
  BCS_out(11:13) = (mu_s * w0_seg1) - w1_seg2
  BCS_out(14)    = NORM2(w1_seg2) - 1.0d0

  !-------------------!
  !     Segment 3     !
  !-------------------!
  ! Boundary conditions - Segments 1 and 3
  BCS_out(15:17) = x1_seg3 - x0_seg1

  !-------------------!
  !     Segment 4     !
  !-------------------!
  BCS_out(18:20) = x0_seg4 - x0_seg3 - (A_perturb * d_perturb)
  BCS_out(21)    = DOT_PRODUCT(x1_seg4 - x0_seg2, w0_seg2)
  ! BCS_out(22)    = NORM2(x1_seg4 - x0_seg2) - eta
  BCS_out(22)    = (NORM2(x1_seg4 - x0_seg2) ** 2.0d0) - eta

END SUBROUTINE BCS_SEGS

SUBROUTINE BCS_SEGS_DFDU(NDIM, NBC, U0, U1, PAR, J_out)
  ! Encoding of the Jacobian of the boundary conditions with
  ! respect to the initial state-space vector components.
  !
  ! Input
  ! ----------
  ! NDIM : INTEGER
  !     Dimension of the algebraic or ODE system
  ! NBC  : INTEGER
  !     Number of boundary conditions
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! PAR : REAL(KIND=8)
  !     Equation parameters
  !
  ! Output
  ! ----------
  ! J_out  : REAL(KIND=8), ARRAY
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
  ! Number of boundary conditions
  INTEGER, INTENT(IN)         :: NBC
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(NDIM)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(NDIM)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: J_out(NBC, 2*NDIM)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
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
  ! Cheaty MATLAB things
  REAL(KIND=8)                :: t2, t3, t4, t5, t6, t7, t8, t9
  REAL(KIND=8)                :: t10, t11, t12, t13, t14, t15
  ! pi
  REAL(KIND=8), PARAMETER     :: pi = 3.14159265358979323846d0

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
  mu_s          = PAR(pdim+5)
  ! Distance from perturbed segment to \Gamma
  eta           = PAR(pdim+6)
  ! Perturbation amplitude
  A_perturb     = PAR(pdim+7)
  ! Angle at which perturbation is applied
  theta_perturb = PAR(pdim+8)
  ! Azimuthal angle at which perturbation is applied
  phi_perturb   = PAR(pdim+9)
  
  !============================================================================!
  !                        CHEATY MATLAB CONSTANT THINGS                       !
  !============================================================================!
  ! Initial vectors
  t2 = ABS(w1_seg2(1))
  t3 = ABS(w1_seg2(2))
  t4 = ABS(w1_seg2(3))
  t5 = x0_seg2(1) * 2.0d0
  t6 = x0_seg2(2) * 2.0d0
  t7 = x0_seg2(3) * 2.0d0
  t8 = x1_seg4(1) * 2.0d0
  t9 = x1_seg4(2) * 2.0d0
  t10 = x1_seg4(3) * 2.0d0
  t11 = t2 ** 2
  t12 = t3 ** 2
  t13 = t4 ** 2
  t14 = t11 + t12 + t13
  t15 = 1.0d0 / SQRT(t14)

  !============================================================================!
  !                             JACOBIAN ENCODING                              !
  !============================================================================!
  !------------------------------------------!
  !     Segments 1 and 2: State - Vector     !
  !------------------------------------------!
  ! dBCS / dU0 - Initial vectors
  J_out(1, 1)   = 1.0d0
  J_out(1, 25)  = -1.0d0
  J_out(2, 2)   = 1.0d0
  J_out(2, 26)  = -1.0d0
  J_out(3, 3)   = 1.0d0
  J_out(3, 27)  = -1.0d0

  J_out(4, 7)   = -1.0d0
  J_out(4, 19)  = 1.0d0
  J_out(5, 8)   = -1.0d0
  J_out(5, 20)  = 1.0d0
  J_out(6, 9)   = -1.0d0
  J_out(6, 21)  = 1.0d0

  J_out(7, 1)   = -gamma * (x0_seg1(3) + 1.0d0)
  J_out(7, 3)   = -gamma * x0_seg1(1)

  J_out(8, 10)  = -1.0d0
  J_out(8, 22)  = 1.0d0
  J_out(9, 11)  = -1.0d0
  J_out(9, 23)  = 1.0d0
  J_out(10, 12) = -1.0d0
  J_out(10, 24) = 1.0d0

  J_out(11, 4)  = mu_s
  J_out(11, 28) = -1.0d0
  J_out(12, 5)  = mu_s
  J_out(12, 29) = -1.0d0
  J_out(13, 6)  = mu_s
  J_out(13, 30) = -1.0d0

  J_out(14, 28) = t2 * t15 * ((w1_seg2(1) / ABS(w1_seg2(1))))
  J_out(14, 29) = t3 * t15 * ((w1_seg2(2) / ABS(w1_seg2(2))))
  J_out(14, 30) = t4 * t15 * ((w1_seg2(3) / ABS(w1_seg2(3))))

  !-------------------!
  !     Segment 3     !
  !-------------------!
  J_out(15, 1)  = -1.0d0
  J_out(15, 31) = 1.0d0
  J_out(16, 2)  = -1.0d0
  J_out(16, 32) = 1.0d0
  J_out(17, 3)  = -1.0d0
  J_out(17, 33) = 1.0d0

  !-------------------!
  !     Segment 4     !
  !-------------------!
  J_out(18, 13) = -1.0d0
  J_out(18, 16) = 1.0d0
  J_out(19, 14) = -1.0d0
  J_out(19, 17) = 1.0d0
  J_out(20, 15) = -1.0d0
  J_out(20, 18) = 1.0d0

  J_out(21, 7)  = -w0_seg2(1)
  J_out(21, 8)  = -w0_seg2(2)
  J_out(21, 9)  = -w0_seg2(3)
  J_out(21, 10) = -x0_seg2(1) + x1_seg4(1)
  J_out(21, 11) = -x0_seg2(2) + x1_seg4(2)
  J_out(21, 12) = -x0_seg2(3) + x1_seg4(3)
  J_out(21, 34) = w0_seg2(1)
  J_out(21, 35) = w0_seg2(2)
  J_out(21, 36) = w0_seg2(3)
  
  J_out(22, 7)  = t5 - t8
  J_out(22, 8)  = t6 - t9
  J_out(22, 9)  = t7 - t10
  J_out(22, 34) = -t5 + t8
  J_out(22, 35) = -t6 + t9
  J_out(22, 36) = -t7 + t10

END SUBROUTINE BCS_SEGS_DFDU

SUBROUTINE BCS_SEGS_DFDP(NDIM, NBC, U0, U1, PAR, J_out)
  ! Encoding of the Jacobian of the boundary conditions with
  ! respect to the parameters.
  !
  ! Input
  ! ----------
  ! NDIM : INTEGER
  !     Dimension of the algebraic or ODE system
  ! NBC  : INTEGER
  !     Number of boundary conditions
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! PAR : REAL(KIND=8)
  !     Equation parameters
  !
  ! Output
  ! ----------
  ! J_out  : REAL(KIND=8), ARRAY
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
  ! Number of boundary conditions
  INTEGER, INTENT(IN)         :: NBC
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(NDIM)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(NDIM)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: J_out(NBC, *)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
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
  ! pi
  REAL(KIND=8), PARAMETER     :: pi = 3.14159265358979323846d0

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
  mu_s          = PAR(pdim+5)
  ! Distance from perturbed segment to \Gamma
  eta           = PAR(pdim+6)
  ! Perturbation amplitude
  A_perturb     = PAR(pdim+7)
  ! Angle at which perturbation is applied
  theta_perturb = PAR(pdim+8)
  ! Azimuthal angle at which perturbation is applied
  phi_perturb   = PAR(pdim+9)

  !============================================================================!
  !                             JACOBIAN ENCODING                              !
  !============================================================================!
  J_out(7, 1)   = A_pump - x0_seg1(1) - (x0_seg1(1) * x0_seg1(3))
  J_out(7, 2)   = gamma
  J_out(11, 9)  = w0_seg1(1)
  J_out(12, 9)  = w0_seg1(2)
  J_out(13, 9)  = w0_seg1(3)
  J_out(18, 11) = -COS(theta_perturb * (2.0d0 * pi))
  J_out(18, 12) = A_perturb * (2.0d0 * pi) * SIN(theta_perturb * (2.0d0 * pi))
  J_out(20, 11) = -SIN(theta_perturb * (2.0d0 * pi))
  J_out(20, 12) = -A_perturb * (2.0d0 * pi) * COS(theta_perturb * (2.0d0 * pi))
  J_out(22, 10) = -1.0d0

END SUBROUTINE BCS_SEGS_DFDP

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