!--------------------------------------------------------------------
! Slightly rewritten
t2 = -x1_seg4(1)
t3 = -x1_seg4(2)
t4 = -x1_seg4(3)
t5 = t2 + x0_seg2(1)
t6 = t3 + x0_seg2(2)
t7 = t4 + x0_seg2(3)

BCS_out(1)  = x0_seg1(1) - x1_seg2(1)
BCS_out(2)  = x0_seg1(2) - x1_seg2(2)
BCS_out(3)  = x0_seg1(3) - x1_seg2(3)
BCS_out(4)  = -x0_seg2(1) + x1_seg1(1)
BCS_out(5)  = -x0_seg2(2) + x1_seg1(2)
BCS_out(6)  = -x0_seg2(3) + x1_seg1(3)
BCS_out(7)  = gamma * (A_pump - x0_seg1(1) - (x0_seg1(1) * x0_seg1(3)))

BCS_out(8)  = w0_seg1(1) - w1_seg2(1)
BCS_out(9)  = w0_seg1(2) - w1_seg2(2)
BCS_out(10) = w0_seg1(3) - w1_seg2(3)
BCS_out(11) = -w1_seg1(1) + mu_s * w0_seg2(1)
BCS_out(12) = -w1_seg1(2) + mu_s * w0_seg2(2)
BCS_out(13) = -w1_seg1(3) + mu_s * w0_seg2(3)
BCS_out(14) = SQRT((w0_seg2(1) ** 2) + (w0_seg2(2) ** 2) + (w0_seg2(3) ** 2)) - 1.0d0

BCS_out(15) = -x0_seg1(1) + x1_seg3(1)
BCS_out(16) = -x0_seg1(2) + x1_seg3(2)
BCS_out(17) = -x0_seg1(3) + x1_seg3(3)

BCS_out(18) = -x0_seg3(1) + x0_seg4(1) - A_perturb * cos(theta_perturb)
BCS_out(19) = -x0_seg3(2) + x0_seg4(2)
BCS_out(20) = -x0_seg3(3) + x0_seg4(3) - A_perturb * sin(theta_perturb)
BCS_out(21) = -t5 * w0_seg2(1) - t6 * w0_seg2(2) - t7 * w0_seg2(3)
BCS_out(22) = -eta + (t5 ** 2) + (t6 ** 2) + (t7 ** 2)
