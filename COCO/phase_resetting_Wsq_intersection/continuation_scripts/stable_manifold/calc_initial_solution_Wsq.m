function data_out = calc_initial_solution_Wsq(run_in, label_in, epsilon_in, funcs_in)
  % data_out = calc_initial_solution_Wsq(run_in, label_in, epsilon_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.
  %
  % Parameters
  % ----------
  % run_in : string
  %     The run identifier for the continuation problem.
  % label_in : int
  %     The solution label for the continuation problem.
  % epsilon_in : float
  %     Initial distance of trajectory segment from the equilibrium point.
  % funcs_in : cell
  %     List of all field function handles.
  %
  % Returns
  % -------
  % data_out : struct
  %     Structure containing the stable manifold data.
  %     Fields:
  %         - x0      : Initial state vector for the stable manifold segment.
  %         - t0      : Initial time.
  %         - ls      : Stable eigenvalue.
  %         - vs      : Stable eigenvector.
  %         - epsilon : Initial distance from the equilibrium.
  %         - p0      : Parameters of the solution.
  %
  % See Also
  % --------
  % ep_read_solution

  %-------------------%
  %     Read Data     %
  %-------------------%
  % x_pos solution
  [sol_xpos, data] = ep_read_solution('xpos',  run_in, label_in);
  xpos     = sol_xpos.x;
  % Parameters
  p        = sol_xpos.p;

  %------------------------------------------------%
  %     Calculate Eigenvectors and Eigenvalues     %
  %------------------------------------------------%
  % Calculate Jacobian of the equilibrium point
  field_DFDX = funcs_in.field{2};
  J = field_DFDX(xpos, p);

  % Calculate eigenvalues and eigenvectors
  [eigval, eigvec] = eig(J);

  % Stable eigenvalue and eigenvector
  ls = eigval(3, 3);
  vs = eigvec(:, 3);

  %-------------------------------------%
  %     Setup Stable Manifold Stuff     %
  %-------------------------------------%
  % Initial time
  t0 = 0;

  % Initial state vector
  x0_Wsq = xpos + (epsilon_in * vs);

  %----------------%
  %     Output     %
  %----------------%
  % Load PO solution data
  data_out.xdim = data.xdim;
  data_out.pdim = data.pdim;
  data_out.x0      = x0_Wsq';
  data_out.t0      = t0;
  data_out.ls      = ls;
  data_out.vs      = vs;
  data_out.epsilon = epsilon_in;
  data_out.p0      = p;

end