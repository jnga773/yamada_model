      t2 = -x1_seg41
      t3 = -x1_seg42
      t4 = -x1_seg43
      t5 = t2+x0_seg21
      t6 = t3+x0_seg22
      t7 = t4+x0_seg23
      A0(1,1) = x0_seg11-x1_seg21
      A0(2,1) = x0_seg12-x1_seg22
      A0(3,1) = x0_seg13-x1_seg23
      A0(4,1) = -x0_seg21+x1_seg11
      A0(5,1) = -x0_seg22+x1_seg12
      A0(6,1) = -x0_seg23+x1_seg13
      A0(7,1) = -gam*(-A_pump+x0_seg11+x0_seg11*x0_seg13)
      A0(8,1) = -w0_seg21+w1_seg11
      A0(9,1) = -w0_seg22+w1_seg12
      A0(10,1) = -w0_seg23+w1_seg13
      A0(11,1) = -w1_seg21+mu_s*w0_seg11
      A0(12,1) = -w1_seg22+mu_s*w0_seg12
      A0(13,1) = -w1_seg23+mu_s*w0_seg13
      A0(14,1) = sqrt(abs(w1_seg21)**2+abs(w1_seg22)**2+abs(w1_seg23)**2
     &)-1.0D0
      A0(15,1) = -x0_seg11+x1_seg31
      A0(16,1) = -x0_seg12+x1_seg32
      A0(17,1) = -x0_seg13+x1_seg33
      A0(18,1) = -x0_seg31+x0_seg41-A_perturb*cos(theta_perturb)
      A0(19,1) = -x0_seg32+x0_seg42
      A0(20,1) = -x0_seg33+x0_seg43-A_perturb*sin(theta_perturb)
      A0(21,1) = -t5*w0_seg21-t6*w0_seg22-t7*w0_seg23
      A0(22,1) = -eta+t5**2+t6**2+t7**2