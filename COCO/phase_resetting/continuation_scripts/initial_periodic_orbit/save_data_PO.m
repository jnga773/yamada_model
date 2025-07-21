function save_data_PO(run_in, label_in, filename_out)
  % save_data_PO(run_in, label_in, filename_out)
  %
  % Reads periodic orbit solution data from COCO solution, calculates the
  % one-dimensional stable manifold of the "central" saddle point 'q', and
  % saves the data to './data/initial_PO.mat'.
  %
  % Parameters
  % ----------
  % run_in : string
  %     The run identifier for the continuation problem.
  % label_in : int
  %     The solution label for the continuation problem.
  % filename_out : string
  %     Filename for the Matlab .mat data file.
  %
  % Returns
  % -------
  % data_out : struct
  %     Structure containing the initial periodic solution data.
  %     Fields:
  %         - xbp_PO : State space solution of the periodic orbit.
  %         - tbp_PO : Temporal data of the periodic orbit.
  %         - T_PO : Period of the periodic orbit.
  %         - x0 : Equilibrium point at x0.
  %         - xpos : Equilibrium point at xpos.
  %         - xneg : Equilibrium point at xneg.
  %         - p : Parameters of the solution.
  %         - pnames : Names of the parameters.
  %         - xdim : Dimension of the state space solution.
  %         - pdim : Dimension of the parameter space.
  %
  % See Also
  % --------
  % coll_read_solution, ep_read_solution

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    label_in double
    filename_out char
  end

  %-----------------------------------%
  %     Read Data: Periodic Orbit     %
  %-----------------------------------%
  % Read COCO solution
  [sol_PO, data_PO] = coll_read_solution('initial_PO', run_in, label_in);

  % State space solution
  xbp_PO = sol_PO.xbp;
  % Time
  tbp_PO = sol_PO.tbp;
  % Period
  T_PO   = sol_PO.T;
  % Parameters
  p      = sol_PO.p;
  pnames = data_PO.pnames;

  %-----------------------------------%
  %     Calculate Stable Manifold     %
  %-----------------------------------%
  Wq_stable = calc_stable_manifold(run_in, label_in);

  %---------------------------------------%
  %     Read Data: Equilibrium Points     %
  %---------------------------------------%
  % COCO solutions
  sol_0   = ep_read_solution('x0', run_in, label_in);
  sol_neg = ep_read_solution('xneg', run_in, label_in);
  sol_pos = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium points
  x0   = sol_0.x;
  xpos = sol_pos.x;
  xneg = sol_neg.x;

  %----------------%
  %     Output     %
  %----------------%
  data_out.xbp_PO  = xbp_PO;
  data_out.tbp_PO  = tbp_PO;
  data_out.T_PO    = T_PO;

  data_out.x0      = x0;
  data_out.xpos    = xpos;
  data_out.xneg    = xneg;

  data_out.p       = p;
  data_out.pnames  = pnames;

  data_out.xdim    = length(x0);
  data_out.pdim    = length(p);

  data_out.Wq_s    = Wq_stable;

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_out, '-struct', 'data_out');

end

function x_out = calc_stable_manifold(run_in, label_in)
  % x_out = calc_stable_manifold(run_in, label_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    label_in double
  end

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read EP solution
  sol_EP   = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium point
  x_pos = sol_EP.x;
  % Parameters
  p     = sol_EP.p;
  
  %------------------------------%
  %     Calculate EigenStuff     %
  %------------------------------%
  % Jacobian
  J_stable = yamada_DFDX(x_pos, p);

  % Calculate eigenvalues and eigenvectors
  [eigvec, eigval] = eig(J_stable);

  % Indices for stable eigenvectors (eigval < 0)
  % stable_index = find(diag(eigval) < 0);
  stable_index = 3;

  % Stable eigenvector
  vec_s = eigvec(:, stable_index);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps1 = -0.01;
  % Time span
  % t_span1 = [0.0, -17.0];
  t_span1 = -16.5 : 0.01 : 0.0;
  t_span1 = flip(t_span1);

  % Initial vector
  x_init1 = x_pos + (eps1 * vec_s);

  % Integrate using ode45
  [~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span1, x_init1);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps2 = 0.01;
  % Time span
  % t_span2 = [0.0, -25.0];
  t_span2 = -25.0 : 0.01 : 0.0;
  t_span2 = flip(t_span2);

  % Initial vector
  x_init2 = x_pos + (eps2 * vec_s);

  % Integrate using ode45
  [~, W2] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span2, x_init2);

  %----------------%
  %     Output     %
  %----------------%
  x_out = [flip(W2); W1];

end