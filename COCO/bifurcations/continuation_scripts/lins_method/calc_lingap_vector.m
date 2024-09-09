function data_out = calc_lingap_vector(run_in, label_in)
  % data_out = calc_lingap_vector(run_in)
  %
  % Calculates the Lin gap vector, the distance between the unstable
  % and stable trajectories, and the Lin phase condition from [run_in].
  %
  % Input
  % ----------
  % run_in : str
  %     The string identifier for the previous COCO run that we will
  %     read information from.
  %
  % Output
  % ----------
  % data_out : structure
  %     Data structure containing the Lin gap vector, distance, and phase.

  %--------------------------------------%
  %     Read Solutions from [run_in]     %
  %--------------------------------------%
  % Extract solution of unstable manifold
  [solu, datau] = coll_read_solution('unstable', run_in, 1);

  % Extract solution of stable manifold
  [sols, datas] = coll_read_solution('stable', run_in, label_in);

  % Final point of unstable manifold
  x1_unstable = solu.xbp(end, :);

  % Initial point of stable manifold
  x0_stable = sols.xbp(1, :);

  %----------------------------------%
  %     Calculate Lin Gap Vector     %
  %----------------------------------%
  % Lin gap vector
  % vgap = x1_unstable - x0_stable;
  vgap = x0_stable - x1_unstable;

  % Calculate norm (i.e., the initial value of lingap)
  vgap_norm = norm(vgap, 2);

  %----------------%
  %     Output     %
  %----------------%
  data_out.xdim    = datau.xdim;
  data_out.pdim    = datau.pdim;
  
  % Add normalised Lin gap vector to data_out
  data_out.vgap    = vgap / vgap_norm;

  % Add value of norm to data_out
  data_out.lingap0 = vgap_norm;
  data_out.lingap  = vgap_norm;

end
