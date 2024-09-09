function bcs_coco_out = bcs_initial_symbolic()
  % bcs_coco_out = bcs_initial_symbolic()
  %
  % COCO compatible encoding for the "initial" boundary conditions of the two
  % trajectory segments. Both segments start near the equilibrium point, x0.
  % The unstable manifold is solved in forwards time, with
  %           x_u(0) = x0 + eps1 * vu.
  %
  % The stable manifold is solved in reverse time, with
  %           x_s(T2) = x0 + eps2(vs1 * cos(theta) + vs2 * sin(theta)) .
  %
  % Here, eps1 and eps2 are the distances from of the trajectories to the equilbrium
  % points, vu is the unstable eigenvector of the Jacobian, vs1 and vs2 are the
  % stable eigenvectors of the Jacobian, and theta is the angle between the two.
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %          u_in(1:3)   - The initial point of the unstable manifold (x0_unstable),
  %          u_in(4:6)   - The final point of the stable manifold (x1_stable),
  %          u_in(7:9)   - The x_neg equilibrium point (x_neg),
  %          u_in(10:12) - Unstable eigenvector,
  %          u_in(13:15) - Stable eigenvector 1,
  %          u_in(16:18) - Stable eigenvector 2,
  %          u_in(19:21) - The epsilon spacings and angle (eps).
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
  % Initial vector of the unstable manifold
  x0_unstable = sym('x0_unstable', [xdim, 1]);
  % Final vector of the stable manifold
  x1_stable   = sym('x1_stable', [xdim, 1]);
  % x_neg equilibrium point
  x_neg       = sym('x_neg', [xdim, 1]);

  % Unstable eigenvector
  vu          = sym('vu', [xdim, 1]);
  % Stable eigenvectors
  vs1         = sym('vs1', [xdim, 1]);
  vs2         = sym('vs2', [xdim, 1]);

  % Epsilon spacings and angle
  syms eps1 eps2 theta

  % Combined vector
  uvec = [x0_unstable;
          x1_stable;
          x_neg;
          vu;
          vs1;
          vs2;
          eps1; eps2; theta];

  %---------------------------------------%
  %     Calculate Boundary Conditions     %
  %---------------------------------------%
  % Unstable boundary condition
  x_init_u = x_neg + (eps1 * vu);

  % Stable boundary condition
  x_final_s = x_neg + eps2 * ((vs1 * cos(theta)) + (vs2 * sin(theta)));

  % Boundary condition equations
  bcs1 = x0_unstable - x_init_u;
  bcs2 = x1_stable - x_final_s;

  % Boundary condition vector
  bcs = [bcs1; bcs2];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_initial';

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