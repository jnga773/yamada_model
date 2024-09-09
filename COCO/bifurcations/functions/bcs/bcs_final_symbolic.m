function bcs_coco_out = bcs_final_symbolic()
  % bcs_coco_out = bcs_final_symbolic()
  %
  % Boundary conditions of the two trajectory segments. Both segments end on
  % the plane \Sigma, which we define by a point and a normal vector in
  % lins_method_setup().
  % The unstable manifold is solved in forwards time, with
  %           x_u(T1) \in \Sigma.
  % The stable manifold is solved in reverse time, with
  %           x_s(0) \in \Sigma.
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %          u_in(1:3) - The final point of the unstable manifold (x1_unstable),
  %          u_in(4:6) - The initial point of the stable manifold (x0_stable),
  %          u_in(7:9) - The x_pos non-trivial equilibrium point in the
  %                      centre of the homoclinic; used to fix the z-position
  %                      of the \Sigma plane (x_pos).
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
  x1_unstable = sym('x1_unstable', [xdim, 1]);
  % Final vector of the stable manifold
  x0_stable   = sym('x0_stable', [xdim, 1]);
  % x_pos equilibrium point
  x_pos       = sym('x_neg', [xdim, 1]);

  % Combined vector
  uvec = [x1_unstable;
          x0_stable;
          x_pos];

  % Normal vector to plane
  normal = [0, 0, 1];
  % Intersection point on plane
  pt0    = [0; 0; x_pos(3)];

  %---------------------------------------%
  %     Calculate Boundary Conditions     %
  %---------------------------------------%
  % Boundary condition equations
  bcs1 = normal * (x1_unstable - pt0);
  bcs2 = normal * (x0_stable - pt0);

  % Boundary condition vector
  bcs = [bcs1; bcs2];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_final';

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