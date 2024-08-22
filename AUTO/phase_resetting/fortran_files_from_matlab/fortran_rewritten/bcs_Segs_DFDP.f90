!--------------------------------------------------------------------
! Slightly rewritten

J_out(7, 1)   = A_pump - x0_seg1(1) - (x0_seg1(1) * x0_seg1(3))
J_out(7, 2)   = gamma
J_out(11, 9)  = w0_seg1(1)
J_out(12, 9)  = w0_seg1(2)
J_out(13, 9)  = w0_seg1(3)
J_out(18, 11) = -COS(theta_perturb)
J_out(18, 12) = A_perturb * SIN(theta_perturb)
J_out(20, 11) = -SIN(theta_perturb)
J_out(20, 12) = -A_perturb * COS(theta_perturb)
J_out(22, 10) = -1.0D0

