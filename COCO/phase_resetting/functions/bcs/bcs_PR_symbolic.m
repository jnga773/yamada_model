function bcs_coco_out = bcs_PR_symbolic()
  % bcs_coco_out = bcs_PR_symbolic()
  %
  % Boundary conditions for the four segments of the phase-resetting problem:
  %                          x1(0) - x2(1) = 0 ,
  %                          x1(1) - x2(0) = 0 ,
  %                          e1 . F(x1(0)) = 0 ,
  %                          w1(0) - w2(1) = 0 ,
  %                     mu_s w2(0) - w1(1) = 0 ,
  %                            |w2(0)| - 1 = 0 ,
  %                          x3(1) - x1(0) = 0 ,
  %                x4(0) - x3(0) - A_p d_p = 0 ,
  %                (x4(1) - x2(0)) . w2(0) = 0 ,
  %               | x4(1) - x2(0) | - \eta = 0 .
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %            u_in(1:3)   - x(0) of segment 1,
  %            u_in(4:6)   - w(0) of segment 1,
  %            u_in(7:9)   - x(0) of segment 2,
  %            u_in(10:12) - w(0) of segment 2,
  %            u_in(13:15) - x(0) of segment 3,
  %            u_in(16:18) - x(0) of segment 4,
  %            u_in(19:21) - x(1) of segment 1,
  %            u_in(22:24) - w(1) of segment 1,
  %            u_in(25:27) - x(1) of segment 2,
  %            u_in(28:30) - w(1) of segment 2,
  %            u_in(31:33) - x(1) of segment 3,
  %            u_in(34:36) - x(1) of segment 4,
  %            u_in(37:50) - Parameters.
  %
  % Output
  % ----------
  % bcs_coco_out : cell of function handles
  %     List of CoCo-ified symbolic functions for the boundary conditions
  %     Jacobian, and Hessian.

  % State-space dimension
  xdim = 3;

  %============================================================================%
  %                              INPUT PARAMETERS                              %
  %============================================================================%
  %--------------------------------%
  %     Input: Initial Vectors     %
  %--------------------------------%
  % Segment 1 - x(0)
  x0_seg1 = sym('x0_seg1', [xdim, 1]);
  % Segment 1 - w(0)
  w0_seg1 = sym('w0_seg1', [xdim, 1]);
  % Segment 2 - x(0)
  x0_seg2 = sym('x0_seg2', [xdim, 1]);
  % Segment 2 - w(0)
  w0_seg2 = sym('w0_seg2', [xdim, 1]);
  % Segment 3 - x(0)
  x0_seg3 = sym('x0_seg3', [xdim, 1]);
  % Segment 4 - x(0)
  x0_seg4 = sym('x0_seg4', [xdim, 1]);

  %------------------------------%
  %     Input: Final Vectors     %
  %------------------------------%
  % Segment 1 - x(1)
  x1_seg1 = sym('x1_seg1', [xdim, 1]);
  % Segment 1 - w(1)
  w1_seg1 = sym('w1_seg1', [xdim, 1]);
  % Segment 2 - x(1)
  x1_seg2 = sym('x1_seg2', [xdim, 1]);
  % Segment 2 - w(1)
  w1_seg2 = sym('w1_seg2', [xdim, 1]);
  % Segment 3 - x(1)
  x1_seg3 = sym('x1_seg3', [xdim, 1]);
  % Segment 4 - x(1)
  x1_seg4 = sym('x1_seg4', [xdim, 1]);

  %---------------------------%
  %     Input: Parameters     %
  %---------------------------%
  % System parameters
  syms gam A B a
  p_sys = [gam; A; B; a];

  % Phase resetting parameters
  syms k theta_old theta_new mu_s
  syms eta A_perturb theta_perturb phi_perturb
  p_PR = [k; theta_old; theta_new; mu_s;
          eta; A_perturb; theta_perturb; phi_perturb];

  %============================================================================%
  %                         BOUNDARY CONDITION ENCODING                        %
  %============================================================================%
  %---------------------------------%
  %     Segment 1 and Segment 2     %
  %---------------------------------%
  % Vector field
  F_vec = yamada_symbolic_field(x0_seg1, p_sys);

  % Boundary Conditions - Segments 1 and 2
  bcs_seg12_1   = x0_seg1 - x1_seg2;
  bcs_seg12_2   = x0_seg2 - x1_seg1;
  bcs_seg12_3   = F_vec(1);

  % Adjoint Boundary Conditions - Segments 1 and 2
  a_bcs_seg12_1 = w1_seg1 - w0_seg2;
  a_bcs_seg12_2 = (mu_s * w0_seg1) - w1_seg2;
  % a_bcs_seg12_3 = norm(w1_seg1) - 1;
  a_bcs_seg12_3 = norm(w1_seg2) - 1;

  %-------------------%
  %     Segment 3     %
  %-------------------%
  % Boundary Conditions - Segment 3
  bcs_seg3 = x1_seg3 - x0_seg1;

  %-------------------%
  %     Segment 4     %
  %-------------------%
  % d_vec = [cos(theta_perturb) * sin(phi_perturb);
  %          sin(theta_perturb) * sin(phi_perturb);
  %          cos(phi_perturb)];
  d_vec = [cos(theta_perturb * (2.0 * pi));
           0.0;
           sin(theta_perturb) * (2.0 * pi)];

  % Boundary Conditions - Segment 4
  bcs_seg4_1 = x0_seg4 - x0_seg3 - (A_perturb * d_vec);
  bcs_seg4_2 = dot(x1_seg4 - x0_seg2, w0_seg2);
  % bcs_seg4_3 = norm(x1_seg4 - x0_seg2) - eta;

  % The last boundary condition has a singularity in the Jacobian for the initial
  % vector, as the norm is zero. We then redfine this boundary condition as the
  % square.
  diff_vec = x1_seg4 - x0_seg2;
  bcs_seg4_3 = (diff_vec(1) ^ 2) + (diff_vec(2) ^ 2) + (diff_vec(3) ^ 2) - eta;
  % bcs_seg4_3 = ((x1_seg4(1) - x0_seg2(1)) ^ 2) + ((x1_seg4(2) - x0_seg2(2)) ^ 2) + ((x1_seg4(3) - x0_seg2(3)) ^ 2) - eta;

  %============================================================================%
  %                                   OUTPUT                                   %
  %============================================================================%
  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Combined vector
  uvec = [x0_seg1; w0_seg1; x0_seg2; w0_seg2; x0_seg3; x0_seg4;
          x1_seg1; w1_seg1; x1_seg2; w1_seg2; x1_seg3; x1_seg4;
          p_sys; p_PR];

  % Boundary conditions vector
  bcs =  [bcs_seg12_1;  bcs_seg12_2; bcs_seg12_3;
          a_bcs_seg12_1; a_bcs_seg12_2; a_bcs_seg12_3;
          bcs_seg3;
          bcs_seg4_1; bcs_seg4_2; bcs_seg4_3];

  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_PR';

  % COCO Function encoding
  bcs_coco = sco_sym2funcs(bcs, {uvec}, {'u'}, 'filename', filename_out);

  % Function to "CoCo-ify" function outputs: [data_in, y_out] = f(prob_in, data_in, u_in)
  cocoify = @(func_in) @(prob_in, data_in, u_in) deal(data_in, func_in(u_in));

  % List of functions
  func_list = {cocoify(bcs_coco('')), cocoify(bcs_coco('u')), cocoify(bcs_coco({'u', 'u'}))};

  %----------------%
  %     Output     %
  %----------------%
  bcs_coco_out = func_list;

end
