!--------------------------------------------------------------------
! Slightly rewritten
t2 = -A_pump
t3 = -B
t4 = theta_old - 1.0d0
t5 = theta_new - 1.0d0
t6 = -x1_vec(1)
t7 = -x2_vec(1)
t8 = T * theta_new * w1_vec(3) * x1_vec(3)
t9 = t6 + x1_vec(2) + 1.0d0
t10 = t7 + x2_vec(2) + 1.0d0
t11 = T * t5 * w2_vec(3) * x2_vec(3)

! Segment 1: State space
A0(1, 1) = -T * gamma * theta_new * (-A_pump + x1_vec(1) + x1_vec(1) * x1_vec(3))
A0(2, 1) = -T * gamma * theta_new * (-B + x1_vec(2) + a * x1_vec(2) * x1_vec(3))
A0(3, 1) = -T * (-x1_vec(1) + x1_vec(2) + 1.0d0) * theta_new * x1_vec(3)
! Segment 1: Adjoint space
A0(4, 1) = -(T * theta_new * w1_vec(3) * x1_vec(3)) + &
         & T * gamma * theta_new * w1_vec(1) * (x1_vec(3) + 1.0d0)
A0(5, 1) = (T * theta_new * w1_vec(3) * x1_vec(3)) + &
         & T * gamma * theta_new * w1_vec(2) * (a * x1_vec(3) + 1.0d0)
A0(6, 1) = T * (-x1_vec(1) + x1_vec(2) + 1.0d0) * theta_new * w1_vec(3) + &
         & T * gamma * theta_new * w1_vec(1) * x1_vec(1) + &
         & T * a * gamma * theta_new * w1_vec(2) * x1_vec(2)

! Segment 2: State space
A0(7, 1) = T * gamma * (theta_new - 1.0d0) * (-A_pump + x2_vec(1) + x2_vec(1) * x2_vec(3))
A0(8, 1) = T * gamma * (theta_new - 1.0d0) * (-B + x2_vec(2) + a * x2_vec(2) * x2_vec(3))
A0(9, 1) = T * (theta_new - 1.0d0) * (-x2_vec(1) + x2_vec(2) + 1.0d0) * x2_vec(3)
! Segment 2: Adjoint space
A0(10, 1) = (T * (theta_new - 1.0d0) * w2_vec(3) * x2_vec(3)) - &
          & T * gamma * (theta_new - 1.0d0) * w2_vec(1) * (x2_vec(3) + 1.0d0)
A0(11, 1) = -(T * (theta_new - 1.0d0) * w2_vec(3) * x2_vec(3)) - &
          & T * gamma * (theta_new - 1.0d0) * w2_vec(2) * (a * x2_vec(3) + 1.0d0)
A0(12, 1) = -T * (theta_new - 1.0d0) * (-x2_vec(1) + x2_vec(2) + 1.0d0) * w2_vec(3) + &
          & T * gamma * (theta_new - 1.0d0) * -x2_vec(1) * w2_vec(1) - &
          & T * a * gamma * (theta_new - 1.0d0) * w2_vec(2) * x2_vec(2)

! Segment 3: State space
A0(13, 1) = T * gamma * (theta_old - 1.0d0) * (-A_pump + x3_vec(1) + x3_vec(1) * x3_vec(3))
A0(14, 1) = T * gamma * (theta_old - 1.0d0) * (-B + x3_vec(2) + a * x3_vec(2) * x3_vec(3))
A0(15, 1) = T * (theta_old - 1.0d0) * x3_vec(3) * (-x3_vec(1) + x3_vec(2) + 1.0d0)

! Segment 4: State space
A0(16, 1) = -T * gamma * k * (-A_pump + x4_vec(1) + x4_vec(1) * x4_vec(3))
A0(17, 1) = -T * gamma * k * (-B + x4_vec(2) + a * x4_vec(2) * x4_vec(3))
A0(18, 1) = -T * k * x4_vec(3) * (-x4_vec(1) + x4_vec(2) + 1.0d0)

!--------------------------------------------------------------------
! Rewritten betterer
! Segment 1: State space
F_out(1)   = gamma * (A_pump - x1_vec(1) - (x1_vec(1) * x1_vec(3)))
F_out(2)   = gamma * (B - x1_vec(2) - (a * x1_vec(2) * x1_vec(3)))
F_out(3)   = x1_vec(3) * (x1_vec(1) - x1_vec(2) - 1.0d0)
F_out(1:3) = F_out(1:3) * T * theta_new

! Segment 1: Adjoint space
F_out(4)   = -(w1_vec(3) * x1_vec(3)) + &
           & gamma * w1_vec(1) * (x1_vec(3) + 1.0d0)
F_out(5)   = (w1_vec(3) * x1_vec(3)) + &
           & gamma * w1_vec(2) * (a * x1_vec(3) + 1.0d0)
F_out(6)   = -(x1_vec(1) - x1_vec(2) - 1.0d0) * w1_vec(3) + &
           & gamma * w1_vec(1) * x1_vec(1) + &
           & a * gamma * w1_vec(2) * x1_vec(2)
F_out(4:6) = F_out(4:6) * T * theta_new

! Segment 2: State space
F_out(7)   = gamma * (A_pump - x2_vec(1) - (x2_vec(1) * x2_vec(3)))
F_out(8)   = gamma * (B - x2_vec(2) - (a * x2_vec(2) * x2_vec(3)))
F_out(9)   = x2_vec(3) * (x2_vec(1) - x2_vec(2) - 1.0d0)
F_out(7:9) = F_out(7:9) * T * (1.0d0 - theta_new)

! Segment 2: Adjoint space
F_out(10)    = -(w2_vec(3) * x2_vec(3)) + &
             & gamma * w2_vec(1) * (x2_vec(3) + 1.0d0)
F_out(11)    = (w2_vec(3) * x2_vec(3)) + &
             & gamma * w2_vec(2) * (a * x2_vec(3) + 1.0d0)
F_out(12)    = (-x2_vec(1) + x2_vec(2) + 1.0d0) * w2_vec(3) - &
             & gamma * (-x2_vec(1)) * w2_vec(1) + &
             & a * gamma * w2_vec(2) * x2_vec(2)
F_out(10:12) = F_out(10:12) * T * (1.0d0 - theta_new)

! Segment 3: State space
F_out(13)    = gamma * (A_pump - x3_vec(1) - (x3_vec(1) * x3_vec(3)))
F_out(14)    = gamma * (B - x3_vec(2) - (a * x3_vec(2) * x3_vec(3)))
F_out(15)    = x3_vec(3) * (x3_vec(1) - x3_vec(2) - 1.0d0)
F_out(13:15) = F_out(13:15) * T * (1.0d0 - theta_old)

! Segment 4: State space
F_out(16)    = gamma * (A_pump - x4_vec(1) - (x4_vec(1) * x4_vec(3)))
F_out(17)    = gamma * (B - x4_vec(2) - (a * x4_vec(2) * x4_vec(3)))
F_out(18)    = x4_vec(3) * (x4_vec(1) - x4_vec(2) - 1.0d0)
F_out(16:18) = F_out(16:18) * T * k
