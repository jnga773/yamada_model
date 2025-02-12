function data_out = calc_initial_solution_Wsq(filename_in)
  % data_out = calc_initial_solution_Wsq(filename_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.
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

  %------------------------------------------------%
  %     Calculate Eigenvectors and Eigenvalues     %
  %------------------------------------------------%
  % Calculate Jacobian of the equilibrium point
  J = yamada_DFDX(data_out.xpos, data_out.p);

  % Calculate eigenvalues and eigenvectors
  [eigval, eigvec] = eig(J);

  % Stable eigenvalue and eigenvector
  ls = eigval(3, 3);
  vs = eigvec(:, 3);

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
  x_init_1 = data_out.xpos + (eps * vs);
  x_init_2 = data_out.xpos - (eps + vs);

  %----------------%
  %     Output     %
  %----------------%
  % Load PO solution data
  data_out.x_init_1 = x_init_1';
  data_out.x_init_2 = x_init_2';
  data_out.t0       = t0;

  data_out.ls       = ls;
  data_out.vs       = vs;

  data_out.eps      = eps;

end