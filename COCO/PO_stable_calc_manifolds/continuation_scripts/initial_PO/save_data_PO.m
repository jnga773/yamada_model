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
  %         - xbp_end : End point of the periodic orbit.
  %         - vec_s : Strong stable eigenvector.
  %         - lam_s : Eigenvalue of the strong stable eigenvector.
  %
  % See Also
  % --------
  % coll_read_solution, ep_read_solution

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

  %----------------------------------%
  %     Read Data: Floquet Stuff     %
  %----------------------------------%
  % End point of periodic orbits
  xbp_end = xbp_PO(end, :);

  % Monodromy matrix
  chart = coco_read_solution('initial_PO.coll.var', run_in, label_in, 'chart');
  data  = coco_read_solution('initial_PO.coll', run_in, label_in, 'data');

  % Create monodrony matrix
  M1 = chart.x(data.coll_var.v1_idx);

  % Get eigenvalues and eigenvectors of the Monodromy matrix
  [floquet_vec, floquet_eig] = eig(M1);

  % Find index for stable eigenvector? < 1
  [lam_s, min_idx] = min(abs(diag(floquet_eig)));

  % Stable eigenvector
  vec_s = floquet_vec(:, min_idx);
  % Stable eigenvalue (Floquet thingie)
  % lam_s = floquet_eig(min_idx, min_idx);

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
  
  data_out.xbp_end = xbp_end;
  data_out.vec_s   = vec_s;
  data_out.lam_s   = lam_s;

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_out, '-struct', 'data_out');

end