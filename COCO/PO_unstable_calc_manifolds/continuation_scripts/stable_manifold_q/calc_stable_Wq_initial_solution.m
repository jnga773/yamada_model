function data_out = calc_stable_Wq_initial_solution(filename_in)
  % data_out = calc_stable_Wq_initial_solution(filename_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read initial solution stuff from filename_in
  load(filename_in);

  %------------------------------------------------%
  %     Calculate Eigenvectors and Eigenvalues     %
  %------------------------------------------------%
  % Calculate non-trivial steady states
  J = yamada_DFDX(xpos, p);

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
  x_init_1 = xpos + (eps * vs);
  x_init_2 = xpos - (eps + vs);

  %----------------%
  %     Output     %
  %----------------%
  % Read periodic orbit solution
  data_out = load(filename_in);

  data_out.xdim     = length(xpos);
  data_out.pdim     = length(pnames);

  data_out.x_init_1 = x_init_1';
  data_out.x_init_2 = x_init_2';
  data_out.t0       = t0;

  data_out.ls       = ls;
  data_out.vs       = vs;

  data_out.eps      = eps;

end