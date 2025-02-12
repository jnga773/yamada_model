function data_out = calc_initial_solution_WsPO(filename_in)
  % data_out = calc_initial_solution_WsPO(filename_in)
  %
  % Calculate the initial solution to the continuation problem, solving for
  % the strong stable manifold of the periodic orbit.
  %
  % Parameters
  % ----------
  % filename_in : string
  %     Matlab .mat data file to read initial solutions from.
  %
  % Returns
  % -------
  % data_out : struct
  %     Structure containing the stable manifold data.
  %     Fields:
  %         - x_init_1 : Initial state vector for the stable manifold (positive direction).
  %         - x_init_2 : Initial state vector for the stable manifold (negative direction).
  %         - t0 : Initial time.
  %         - ls : Stable eigenvalue.
  %         - vs : Stable eigenvector.
  %         - eps : Initial distance from the equilibrium.
  %         - p : Parameters of the solution.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Load data
  data_out = load(filename_in);

  %-------------------------------------%
  %     Setup Stable Manifold Stuff     %
  %-------------------------------------%
  % Initial distances from the equilibria, along the tangent spaces of the
  % unstable and stable manifolds, to the initial points along the corresponding
  % trajectory segments.
  eps = 0.001;

  % Initial time
  t0 = 0;

  % Initial state vector
  x_init_1 = data_out.xbp_end' + (eps * data_out.vec_s);
  x_init_2 = data_out.xbp_end' - (eps * data_out.vec_s);

  %----------------%
  %     Output     %
  %----------------%
  % Load PO solution data
  data_out.x_init_1 = x_init_1';
  data_out.x_init_2 = x_init_2';
  data_out.t0       = t0;

  data_out.eps      = eps;

end