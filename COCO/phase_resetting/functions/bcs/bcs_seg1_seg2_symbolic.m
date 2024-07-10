function bcs_coco_out = bcs_seg1_seg2_symbolic()
  % bcs_coco_out = bcs_seg1_seg2_symbolic()
  %
  % Boundary conditions for segments 1 and 2 of the phase reset
  % segments:
  %                        x1(0) - x2(1) = 0 ,
  %                        x1(1) - x2(0) = 0 ,
  %                        e1 . F(x1(0)) = 0 ,
  % and the adjoint boundary conditions:
  %                        x1(0) - x2(1) = 0 ,
  %                        x1(1) - x2(0) = 0 ,
  %                        e1 . F(x1(0)) = 0 .
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %           u_in(1:3)   - x(0) of segment 1,
  %           u_in(4:6)   - w(0) of segment 1,
  %           u_in(7:9)   - x(0) of segment 2,
  %           u_in(10:12) - w(0) of segment 2,
  %           u_in(13:15) - x(1) of segment 1,
  %           u_in(16:18) - w(1) of segment 1,
  %           u_in(19:21) - x(1) of segment 2,
  %           u_in(22:24) - w(1) of segment 2,
  %           u_in(25:38) - Parameters.
  %
  % Output
  % ----------
  % bcs_coco_out : cell of function handles
  %     List of CoCo-ified symbolic functions for the boundary conditions
  %     Jacobian, and Hessian.

  % State-space dimension
  xdim = 3;

  %---------------%
  %     Input     %
  %---------------%
  % Segment 1 - x(0)
  x0_seg1 = sym('x0_seg1', [xdim, 1]);

  % Segment 1 - w(0)
  w0_seg1 = sym('w0_seg1', [xdim, 1]);

  % Segment 2 - x(0)
  x0_seg2 = sym('x0_seg2', [xdim, 1]);

  % Segment 2 - w(0)
  w0_seg2 = sym('w0_seg2', [xdim, 1]);

  % Segment 1 - x(1)
  x1_seg1 = sym('x1_seg1', [xdim, 1]);

  % Segment 1 - w(1)
  w1_seg1 = sym('w1_seg1', [xdim, 1]);

  % Segment 2 - x(1)
  x1_seg2 = sym('x1_seg2', [xdim, 1]);

  % Segment 2 - w(1)
  w1_seg2 = sym('w1_seg2', [xdim, 1]);

  % System parameters
  syms gam A B a
  p_sys = [gam; A; B; a];

  % Phase resetting parameters
  syms T k mu_s eta
  syms theta_old theta_new
  syms theta_perturb phi_perturb A_perturb
  p_PR = [T; k; mu_s; eta;
          theta_old; theta_new;
          theta_perturb; phi_perturb; A_perturb];

  % Combined vector
  uvec = [x0_seg1; w0_seg1;
          x0_seg2; w0_seg2;
          x1_seg1; w1_seg1;
          x1_seg2; w1_seg2;
          p_sys; p_PR];

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector field
  F_vec = yamada_symbolic_field(x0_seg1, p_sys);

  % Boundary Conditions - Segments 1 and 2
  bcs_seg12_1   = x0_seg1 - x1_seg2;
  bcs_seg12_2   = x1_seg1 - x0_seg2;
  bcs_seg12_3   = F_vec(1);

  % Adjoint Boundary Conditions - Segments 1 and 2
  a_bcs_seg12_1 = w0_seg1 - w1_seg2;
  a_bcs_seg12_2 = (mu_s * w0_seg2) - w1_seg1;
  a_bcs_seg12_3 = norm(w0_seg2) - 1;

  % Boundary condition vector
  bcs = [bcs_seg12_1; bcs_seg12_2; bcs_seg12_3;
         a_bcs_seg12_1; a_bcs_seg12_2; a_bcs_seg12_3];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_seg1_seg2';

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
