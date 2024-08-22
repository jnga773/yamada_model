!--------------------------------------------------------------------
! Slightly rewritten
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
