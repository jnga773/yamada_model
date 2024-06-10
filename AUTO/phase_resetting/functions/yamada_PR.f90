!==============================================================================!
!                           AUTO-07: PHASE RESETTING                           !
!==============================================================================!
! This module file contains the subroutines neeed to calculate the phase
! resetting problem of the Yamada model in AUTO.
!
! This module contains the following subroutines:
! - FIELD      : Vector field encoding of the Yamada model.
! - FIELD_DFDX : Encoding of the state-space Jacobian.
! - FIELD_DFDP : Encoding of the parameter-space Jacobian.
!
! - FUNC           : The encoded vector field of the continuation problem.
! - FUNC_SEG1      : Vector field encoding of the variational problem of
!                    'Segment 1' of the phase resetting problem.
! - FUNC_SEG2      : Vector field encoding of the variational problem of
!                    'Segment 2' of the phase resetting problem.
! - FUNC_SEG3      : Vector field encoding of the variational problem of
!                    'Segment 3' of the phase resetting problem.
! - FUNC_SEG4      : Vector field encoding of the variational problem of
!                    'Segment 4' of the phase resetting problem.
!
! - BCND           : Define the various boundary conditions (if necessary).
! - BCS_SEG1_SEG2  : 'Segment 1' and 'Segment 2': Boundary conditions for the
!                    phase resetting problem.
! - BCS_SEG3       : 'Segment 3': Boundary conditions for the phase resetting
!                    problem.
! - BCS_SEG4       : 'Segment 4': Boundary conditions for the phase resetting
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
  ! Segment vectors
  REAL(KIND=8)                :: U_seg1(2*xdim)
  REAL(KIND=8)                :: U_seg2(2*xdim)
  REAL(KIND=8)                :: U_seg3(xdim)
  REAL(KIND=8)                :: U_seg4(xdim)

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! Segment 1
  U_seg1 = U(1 : 2*xdim)
  ! Segment 2
  U_seg2 = U(2*xdim+1 : 4*xdim)
  ! Segment 3
  U_seg3 = U(4*xdim+1 : 5*xdim)
  ! Segment 4
  U_seg4 = U(5*xdim+1 : 6*xdim)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  !-------------------!
  !     Segment 1     !
  !-------------------!
  CALL FUNC_SEG1(U_seg1, PAR, F_out(1 : 2*xdim))

  !-------------------!
  !     Segment 2     !
  !-------------------!
  CALL FUNC_SEG2(U_seg2, PAR, F_out(2*xdim+1 : 4*xdim))

  !-------------------!
  !     Segment 3     !
  !-------------------!
  CALL FUNC_SEG3(U_seg3, PAR, F_out(4*xdim+1 : 5*xdim))

  !-------------------!
  !     Segment 4     !
  !-------------------!
  CALL FUNC_SEG4(U_seg4, PAR, F_out(5*xdim+1 : 6*xdim))

  ! End of subroutine
END SUBROUTINE FUNC

SUBROUTINE FUNC_SEG1(U, PAR, F_out)
  ! Vector field encoding for the vector field corresponding to
  ! segment 1 of the phase resetting curve.
  ! Segment 1 goes from gamma_0 to theta_new.
  !
  ! Input
  ! ----------
  ! xdim   : INTEGER
  !     Original dimension of the vector field.
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
  ! Import global variables
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(*)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: F_out(2*xdim)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Original vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! State space Jacobian
  REAL(KIND=8)                :: J_DFDX(xdim, xdim)
  ! State space variables
  REAL(KIND=8)                :: x_vec(xdim)
  ! Perpendicular vectors
  REAL(KIND=8)                :: w_vec(xdim)
  ! Parameters
  REAL(KIND=8)                :: T, k, theta_new

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! State space variables
  x_vec     = U(1 : xdim)
  ! Perpeindicular vectors
  w_vec     = U(xdim+1 : 2*xdim)

  ! Period of orbit
  T         = PAR(pdim+1)
  ! Integer number of perioids
  ! k         = PAR(pdim+2)
  k         = 1.0d0
  ! Phase where segments comes back to \Gamma
  theta_new = PAR(pdim+6)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  ! Calculate vector field
  CALL FIELD(x_vec, PAR, vec_field)

  ! Calculate Jacobian
  CALL FIELD_DFDX(x_vec, PAR, J_DFDX)

  ! Vector field
  F_out(1 : xdim) = k * T * theta_new * vec_field

  ! Adjoint equation
  F_out(xdim+1 : 2*xdim) = -k * T * theta_new * MATMUL(TRANSPOSE(J_DFDX), w_vec)

  ! End of subroutine
END SUBROUTINE FUNC_SEG1

SUBROUTINE FUNC_SEG2(U, PAR, F_out)
  ! Vector field encoding for the vector field corresponding to
  ! segment 2 of the phase resetting curve.
  ! Segment 2 goes from theta_new to gamma_0.
  !
  ! Input
  ! ----------
  ! xdim   : INTEGER
  !     Original dimension of the vector field.
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
  ! Import global variables
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(2*xdim)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: F_out(2*xdim)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Original vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! State space Jacobian
  REAL(KIND=8)                :: J_DFDX(xdim, xdim)
  ! State space variables
  REAL(KIND=8)                :: x_vec(xdim)
  ! Perpendicular vectors
  REAL(KIND=8)                :: w_vec(xdim)
  ! Parameters
  REAL(KIND=8)                :: T, k, theta_new

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! State space variables
  x_vec     = U(1 : xdim)
  ! Perpeindicular vectors
  w_vec     = U(xdim+1 : 2*xdim)

  ! Period of orbit
  T         = PAR(pdim+1)
  ! Integer number of perioids
  ! k         = PAR(pdim+2)
  k         = 1.0d0
  ! Phase where segments comes back to \Gamma
  theta_new = PAR(pdim+6)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  ! Calculate vector field
  CALL FIELD(x_vec, PAR, vec_field)

  ! Calculate Jacobian
  CALL FIELD_DFDX(x_vec, PAR, J_DFDX)

  ! Vector field
  F_out(1 : xdim) = k * T * (1.0d0 - theta_new) * vec_field

  ! Adjoint equation
  F_out(xdim+1 : 2*xdim) = -k * T * (1.0d0 - theta_new) * MATMUL(TRANSPOSE(J_DFDX), w_vec)

  ! End of subroutine
END SUBROUTINE FUNC_SEG2

SUBROUTINE FUNC_SEG3(U, PAR, F_out)
  ! Vector field encoding for the vector field corresponding to
  ! segment 3 of the phase resetting curve.
  ! Segment 3 goes from gamma_0 to theta_old.
  !
  ! Input
  ! ----------
  ! xdim   : INTEGER
  !     Original dimension of the vector field.
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
  ! Import global variables
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(xdim)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: F_out(xdim)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Original vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! State space variables
  REAL(KIND=8)                :: x_vec(xdim)
  ! Parameters
  REAL(KIND=8)                :: T, theta_old

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! State space variables
  x_vec     = U(1:xdim)

  ! Period of orbit
  T         = PAR(pdim+1)
  ! Phase where the perturbation is applied
  theta_old = PAR(pdim+5)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  ! Calculate vector field
  CALL FIELD(x_vec, PAR, vec_field)

  ! Vector field
  F_out(1:xdim) = T * (1.0d0 - theta_old) * vec_field

  ! End of subroutine
END SUBROUTINE FUNC_SEG3

SUBROUTINE FUNC_SEG4(U, PAR, F_out)
  ! Vector field encoding for the vector field corresponding to
  ! segment 4 of the phase resetting curve.
  ! Segment 4 goes from theta_old to theta_new.
  !
  ! Input
  ! ----------
  ! xdim   : INTEGER
  !     Original dimension of the vector field.
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
  ! Import global variables
  USE VECTOR_FIELD

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! State variables of the vector field
  REAL(KIND=8), INTENT(IN)    :: U(xdim)
  ! Equation parameters
  REAL(KIND=8), INTENT(IN)    :: PAR(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Equation or ODE right-hand-side values
  REAL(KIND=8), INTENT(OUT)   :: F_out(xdim)

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Original vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! State space variables
  REAL(KIND=8)                :: x_vec(xdim)
  ! Parameters
  REAL(KIND=8)                :: T, k

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! State space variables
  x_vec     = U(1:xdim)

  ! Period of orbit
  T         = PAR(pdim+1)
  ! Integer number of perioids
  k         = PAR(pdim+2)

  !============================================================================!
  !                            VECTOR FIELD ENCODING                           !
  !============================================================================!
  ! Calculate vector field
  CALL FIELD(x_vec, PAR, vec_field)

  ! Vector field
  F_out(1:xdim) = k * T * vec_field

  ! End of subroutine
END SUBROUTINE FUNC_SEG4

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
  ! Import global variables
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
  !---------------------------------!
  !     Segment 1 and Segment 2     !
  !---------------------------------!
  CALL BCS_SEG1_SEG2(U0, U1, PAR, FB(1 : 4*xdim+2))

  !-------------------!
  !     Segment 3     !
  !-------------------!
  CALL BCS_SEG3(U0, U1, PAR, FB(4*xdim+3 : 5*xdim+2))

  !-------------------!
  !     Segment 4     !
  !-------------------!
  CALL BCS_SEG4(U0, U1, PAR, FB(5*xdim+3 : 6*xdim+4))

END SUBROUTINE BCND

SUBROUTINE BCS_SEG1_SEG2(U0, U1, PAR, FB)
  ! Boundary conditions for segments 1 and 2 of the phase reset
  ! segments:
  !                        x1(0) - x2(1) = 0 ,
  !                        x1(1) - x2(0) = 0 ,
  !                        e1 . F(x1(0)) = 0 ,
  ! and the adjoint boundary conditions:
  !                        w1(0) - w2(1) = 0 ,
  !                   mu_s w1(1) - w2(0) = 0 ,
  !                          |w2(0)| - 1 = 0 ,
  !
  ! Input
  ! ----------
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! PAR  : REAL(KIND=8), ARRAY
  !     Array of the equation parameters
  !
  ! Output
  ! ----------
  ! FB   : REAL(KIND=8), ARRAY
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
  ! Equation parameters.
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(*)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the boundary condition functions
  REAL(KIND=8), INTENT(OUT)   :: FB(4*xdim+2)

  !---------------------------!
  !     SUBROUTINE THINGS     !
  !---------------------------!
  ! The vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! Segment 1: Initial and final vectors
  REAL(KIND=8)                :: x0_seg1(xdim), w0_seg1(xdim)
  REAL(KIND=8)                :: x1_seg1(xdim), w1_seg1(xdim)
  ! Segment 2: Initial and final vectors
  REAL(KIND=8)                :: x0_seg2(xdim), w0_seg2(xdim)
  REAL(KIND=8)                :: x1_seg2(xdim), w1_seg2(xdim)
  ! Parameters
  REAL(KIND=8)                :: mu_s

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! Segment 1: Initial vectors
  x0_seg1 = U0(1 : xdim)
  w0_seg1 = U0(xdim+1 : 2*xdim)

  ! Segment 2: Initial vectors
  x0_seg2 = U0(2*xdim+1 : 3*xdim)
  w0_seg2 = U0(3*xdim+1 : 4*xdim)

  ! Segment 1: Final vectors
  x1_seg1 = U1(1 : xdim)
  w1_seg1 = U1(xdim+1 : 2*xdim)

  ! Segment 2: Final vectors
  x1_seg2 = U1(2*xdim+1 : 3*xdim)
  w1_seg2 = U1(3*xdim+1 : 4*xdim)

  ! Stable Floquet eigenvalue
  mu_s = PAR(pdim+3)

  ! Calculate the vector field
  CALL FIELD(x0_seg1, PAR, vec_field)

  !============================================================================!
  !                        BOUNDARY CONDITION ENCODING                         !
  !============================================================================!
  ! Boundary conditions - Segments 1 and 2
  FB(1 : xdim)            = x0_seg1 - x1_seg2
  FB(xdim+1 : 2*xdim)     = x1_seg1 - x0_seg2
  FB(2*xdim+1)            = vec_field(1)

  ! Adjoint Boundary Conditions - Segments 1 and 2
  FB(2*xdim+2 : 3*xdim+1) = w0_seg1 - w1_seg2
  FB(3*xdim+2 : 4*xdim+1) = (mu_s * w0_seg2) - w1_seg1
  FB(4*xdim+2)            = NORM2(w0_seg2) - 1.0d0

END SUBROUTINE BCS_SEG1_SEG2

SUBROUTINE BCS_SEG3(U0, U1, PAR, FB)
  ! Boundary conditions for segment three of the phase reset
  ! segments:
  !                        x3(1) - x1(0) = 0 .
  !
  ! Input
  ! ----------
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! PAR  : REAL(KIND=8), ARRAY
  !     Array of the equation parameters
  !
  ! Output
  ! ----------
  ! FB   : REAL(KIND=8), ARRAY
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
  ! Equation parameters.
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(*)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the boundary condition functions
  REAL(KIND=8), INTENT(OUT)   :: FB(xdim)

  !---------------------------!
  !     SUBROUTINE THINGS     !
  !---------------------------!
  ! The vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! Segment 1: Initial vector
  REAL(KIND=8)                :: x0_seg1(xdim)
  ! Segment 3: Final vectors
  REAL(KIND=8)                :: x1_seg3(xdim)

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! Segment 1: Initial vectors
  x0_seg1 = U0(1 : xdim)

  ! Segment 2: Final vectors
  x1_seg3 = U1(4*xdim+1 : 5*xdim)

  !============================================================================!
  !                        BOUNDARY CONDITION ENCODING                         !
  !============================================================================!
  ! Boundary conditions - Segments 1 and 3
  FB(1 : xdim) = x1_seg3 - x0_seg1

END SUBROUTINE BCS_SEG3

SUBROUTINE BCS_SEG4(U0, U1, PAR, FB)
  ! Boundary conditions for segment four of the phase reset
  ! segments:
  !                x4(0) - x3(0) - A d_r = 0 ,
  !              (x4(1) - x2(0)) . w2(0) = 0 ,
  !             | x4(1) - x2(0) | - \eta = 0 .
  !
  ! Input
  ! ----------
  ! U0   : REAL(KIND=8), ARRAY
  !     State variable values at the 'left' boundary
  ! U1   : REAL(KIND=8), ARRAY
  !     State variable values at the 'right' boundary
  ! PAR  : REAL(KIND=8), ARRAY
  !     Array of the equation parameters
  !
  ! Output
  ! ----------
  ! FB   : REAL(KIND=8), ARRAY
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
  ! Equation parameters.
  REAL(KIND=8), INTENT(IN)    :: PAR(*)
  ! State variable values at the 'left' boundary
  REAL(KIND=8), INTENT(IN)    :: U0(*)
  ! State variable values at the 'right' boundary
  REAL(KIND=8), INTENT(IN)    :: U1(*)

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! The value of the boundary condition functions
  REAL(KIND=8), INTENT(OUT)   :: FB(xdim+2)

  !---------------------------!
  !     SUBROUTINE THINGS     !
  !---------------------------!
  ! The vector field
  REAL(KIND=8)                :: vec_field(xdim)
  ! Segment 2: Initial vectors
  REAL(KIND=8)                :: x0_seg2(xdim), w0_seg2(xdim)
  ! Segment 3: Initial vector
  REAL(KIND=8)                :: x0_seg3(xdim)
  ! Segmenty 4: Initial and final vectors
  REAL(KIND=8)                :: x0_seg4(xdim), x1_seg4(xdim)
  ! Parameters
  REAL(KIND=8)                :: eta, A_perturb, theta_perturb, phi_perturb
  ! Directional vector
  REAL(KIND=8)                :: d_perturb(xdim)

  !============================================================================!
  !                              INPUT PARAMETERS                              !
  !============================================================================!
  ! Segment 2 - x(0)
  x0_seg2 = U0(2*xdim+1 : 3*xdim)
  ! Segment 2 - w(0)
  w0_seg2 = U0(3*xdim+1 : 4*xdim)

  ! Segment 3 - x(0)
  x0_seg3 = U0(4*xdim+1 : 5*xdim)

  ! Segment 4 - x(0)
  x0_seg4 = U0(5*xdim+1 : 6*xdim)
  ! Segment 4 - x(1)
  x1_seg4 = U1(5*xdim+1 : 6*xdim)

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
  ! Boundary conditions - Segments 1 and 3
  FB(1 : xdim) = x0_seg4 - x0_seg3 - (A_perturb * d_perturb)
  FB(xdim+1)   = DOT_PRODUCT(x1_seg4 - x0_seg2, w0_seg2)
  FB(xdim+2)   = NORM2(x1_seg4 - x0_seg2) - eta

END SUBROUTINE BCS_SEG4

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