function bcs_coco_out = bcs_seg4_symbolic()
  % bcs_coco_out = bcs_seg4_symbolic()
  %
  % Boundary conditions for segment four of the phase reset
  % segments:
  %                x4(0) - x3(0) - A d_r = 0 ,
  %              (x4(1) - x2(0)) . w2(0) = 0 ,
  %             | x4(1) - x2(0) | - \eta = 0 .
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %           u_in(1:3)   - x(0) of segment 2,
  %           u_in(4:6)   - w(0) of segment 2,
  %           u_in(7:9)   - x(0) of segment 3,
  %           u_in(10:12) - x(0) of segment 4,
  %           u_in(13:15) - x(1) of segment 4,
  %           u_in(16:29) - Parameters.
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
  % Segment 2 - x(0)
  x0_seg2 = sym('x0_seg2', [xdim, 1]);

  % Segment 2 - w(0)
  w0_seg2 = sym('w0_seg2', [xdim, 1]);

  % Segment 3 - x(0)
  x0_seg3 = sym('x0_seg3', [xdim, 1]);

  % Segment 4 - x(0)
  x0_seg4 = sym('x0_seg4', [xdim, 1]);

  % Segment 4 - x(1)
  x1_seg4 = sym('x1_seg4', [xdim, 1]);

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
  uvec = [x0_seg2; w0_seg2;
          x0_seg3;
          x0_seg4; x1_seg4;
          p_sys; p_PR];

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Displacement vector
%   d_vec = [sin(phi_perturb) * cos(theta_perturb);
%            sin(phi_perturb) * sin(theta_perturb);
%            cos(phi_perturb)];
  d_vec = [cos(theta_perturb);
           0.0;
           sin(theta_perturb)];

  % Boundary Conditions - Segment 4
  bcs_seg4_1 = x0_seg4 - x0_seg3 - (A_perturb * d_vec);
  bcs_seg4_2 = dot(x1_seg4 - x0_seg2, w0_seg2);
  % bcs_seg4_3 = norm(x1_seg4 - x0_seg2) - eta;
  % bcs_seg4_3 = (norm(x1_seg4 - x0_seg2) ^ 2) - eta;
  bcs_seg4_3 = ((x1_seg4(1) - x0_seg2(1)) ^ 2) + ((x1_seg4(2) - x0_seg2(2)) ^ 2) + ((x1_seg4(3) - x0_seg2(3)) ^ 2) - eta;

  % Boundary condition vector
  bcs = [bcs_seg4_1;
         bcs_seg4_2;
         bcs_seg4_3];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_seg4';

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
