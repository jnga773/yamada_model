function save_initial_PO_data(run_in, label_in)
  % save_initial_PO_data(run_in, label_in)
  %
  % Reads periodic orbit solution data from COCO solution, calculates the
  % one-dimensional stable manifold of the "central" saddle point 'q', and
  % saves the data to './data/initial_PO.mat'.

  % Data matrix filename
  filename_out = './data_mat/initial_PO.mat';

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Periodic orbit solution
  sol_PO = coll_read_solution('initial_PO', run_in, label_in);
  xbp_PO = sol_PO.xbp;
  p_PO   = sol_PO.p;

  % Equilibrium point solutions
  sol_pos = ep_read_solution('xpos', run_in, label_in);
  sol_neg = ep_read_solution('xneg', run_in, label_in);
  sol_0   = ep_read_solution('x0', run_in, label_in);

  x_pos = sol_pos.x;
  x_neg = sol_neg.x;
  x_0   = sol_0.x;

  tbp_PO = sol_PO.tbp;

  %-----------------------------------%
  %     Calculate Stable Manifold     %
  %-----------------------------------%
  W_q_stable = calc_stable_manifold(run_in, label_in);

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_out, 'x_pos', 'x_neg', 'x_0', 'xbp_PO', 'p_PO', 'W_q_stable', 'tbp_PO');

end