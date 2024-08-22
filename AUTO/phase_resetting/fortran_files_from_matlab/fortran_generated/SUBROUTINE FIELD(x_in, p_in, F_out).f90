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

  END SUBROUTINE ADJOINT

  SUBROUTINE ADJOINT_DFDX(x_in, p_in, J_out)
    ! Encoding of the state-space Jacobian of the adjoint equations.
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
    ! Grab the state-space variables from x_in
    G  = x_in(1)
    Q  = x_in(2)
    I  = x_in(3)
    ! Grab the adjoint-space variables from x_in
    w1 = x_in(4)
    w2 = x_in(5)
    w3 = x_in(6)

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

  SUBROUTINE ADJOINT_DFDP(x_in, p_in, J_out)
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

  END SUBROUTINE ADJOINT_DFDP