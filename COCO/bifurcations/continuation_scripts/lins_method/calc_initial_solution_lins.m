function data_out = calc_initial_solution_lins(run_in, label_in)
  % data_out = calc_initial_solution_lins(run_in, label_in)
  %
  % Reads data from previous run solution and calculates the 
  % initial conditions for the various different trajectory segments.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read COCO solution for periodic orbit
  [sol_PO, data_PO] = coll_read_solution('homo.po.orb', run_in, label_in);

  % Parameters
  p = sol_PO.p;
  % Solution
  xbp_approx = sol_PO.xbp;
  tbp_approx = sol_PO.tbp;
  T_approx   = sol_PO.T;

  % Parameter names
  pnames = data_PO.pnames;

  % Calculate other equilibrium point
  [xpos, xneg] = non_trivial_ss(p);

  % Dimensions
  xdim = length(xneg);
  pdim = length(p);

  %------------------------------------------------%
  %     Calculate Eigenvectors and Eigenvalues     %
  %------------------------------------------------%
  % Calculate non-trivial steady states
  [vu, vs1, vs2, eigvals] = unstable_stable_eigenvectors(xneg, p);

  % Eigenvectors
  lu  = eigvals(1);
  ls1 = eigvals(2);
  ls2 = eigvals(3);

  %----------------------------------%
  %     Setup Lin's Method Stuff     %
  %----------------------------------%
  % Initial distances from the equilibria, along the tangent spaces of the
  % unstable and stable manifolds, to the initial points along the corresponding
  % trajectory segments.
  eps1 = 0.009;
  eps2 = 0.005;

  % Angle for stable vector component
  theta0 = 5.969;

  % Lin epsilons vector
  epsilon0 = [eps1; eps2; theta0];

  %-----------------------------%
  %     Boundary Conditions     %
  %-----------------------------%
  % Normal vector to hyperplane \Sigma (just the y-axis at x=0.5)
  normal = [0, 0, 1];

  % Intersection point for hyperplane
  pt0 = [1.25; 0.66; xpos(3)];

  % Initial time
  t0 = 0;

  % Unstable Manifold: Initial point
  x_init_u = xneg' + (eps1 * vu');
  % Unstable Manifold: Final point
  x_final_u = pt0;

  % Stable Manifold: Initial point
  x_init_s = xneg' + eps2 * (cos(theta0) * vs1' + sin(theta0) * vs2');
  % Stable Manifold: Final point
  x_final_s = pt0;

  %----------------%
  %     Output     %
  %----------------%
  data_out.xdim         = xdim;
  data_out.pdim         = pdim;
  data_out.p0           = p;
  data_out.pnames       = pnames;

  data_out.run_approx   = run_in;
  data_out.label_approx = label_in;
  data_out.xbp_approx   = xbp_approx;
  data_out.tbp_approx   = tbp_approx;
  data_out.T_approx     = T_approx; 
  
  data_out.xpos         = xpos;
  data_out.xneg         = xneg;

  data_out.vu           = vu;
  data_out.lu           = lu;
  data_out.vs1          = vs1;
  data_out.ls1          = ls1;
  data_out.vs2          = vs2;
  data_out.ls2          = ls2;
  
  data_out.normal       = normal;
  data_out.pt0          = pt0;
  data_out.epsilon      = epsilon0;

  data_out.x_init_u     = x_init_u;
  data_out.x_final_u    = x_final_u;

  data_out.t0           = t0;
  data_out.x_init_s     = x_init_s;
  data_out.x_final_s    = x_final_s;

end