function save_stable_W_PO_data(run_in, filename_out)
  % save_data_Wsq(run_in, label_in, filename_out)
  %
  % Reads periodic orbit solution data from COCO solution, and the stable
  % manifold of the stationary point q, and saves to [filename_out].
  %
  % Parameters
  % ----------
  % run_in : string
  %     The run identifier for the continuation problem.
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
  %         - xbp_end : End point of the periodic orbit.
  %         - vec_s : Strong stable eigenvector.
  %         - lam_s : Eigenvalue of the strong stable eigenvector.
  %         - WsPO_1 : One branch of the strong stable manifold.
  %         - WsPO_2 : The other branch of the strong stable manifold.
  %
  % See Also
  % --------
  % coll_read_solution, ep_read_solution

  %-----------------------------------%
  %     Read Data: Periodic Orbit     %
  %-----------------------------------%
  % Read COCO solution
  [sol_PO, data_PO] = coll_read_solution('initial_PO', run_in, 1);

  % State space solution
  xbp_PO = sol_PO.xbp;
  % Time
  tbp_PO = sol_PO.tbp;
  % Period
  T_PO   = sol_PO.T;
  % Parameters
  p      = sol_PO.p;
  pnames = data_PO.pnames;

  %---------------------------------------%
  %     Read Data: Equilibrium Points     %
  %---------------------------------------%
  % COCO solutions
  sol_0   = ep_read_solution('x0', run_in, 1);
  sol_neg = ep_read_solution('xneg', run_in, 1);
  sol_pos = ep_read_solution('xpos', run_in, 1);

  % Equilibrium points
  x0   = sol_0.x;
  xpos = sol_pos.x;
  xneg = sol_neg.x;

  %------------------------------------%
  %     Read Data: Stable Manifold     %
  %------------------------------------%
  % Create empty data arrays
  W1_iso1 = []; W2_iso1 = [];
  W1_iso2 = []; W2_iso2 = [];
  W1_iso3 = []; W2_iso3 = [];

  % Cycle through stable manifold solutions
  bd = coco_bd_read(run_in);
  for label = coco_bd_labs(bd)
    % Grab solution
    [sol1, ~] = coll_read_solution('W1', run_in, label);
    [sol2, ~] = coll_read_solution('W2', run_in, label);

    % Append to data arrays
    W1_iso1 = [W1_iso1, sol1.xbp(:, 1)];
    W1_iso2 = [W1_iso2, sol1.xbp(:, 2)];
    W1_iso3 = [W1_iso3, sol1.xbp(:, 3)];

    W2_iso1 = [W2_iso1, sol2.xbp(:, 1)];
    W2_iso2 = [W2_iso2, sol2.xbp(:, 2)];
    W2_iso3 = [W2_iso3, sol2.xbp(:, 3)];
  end

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

  data_out.WsPO_1 = {W1_iso1, W1_iso2, W1_iso3};
  data_out.WsPO_2 = {W2_iso1, W2_iso2, W2_iso3};

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_out, '-struct', 'data_out');

end