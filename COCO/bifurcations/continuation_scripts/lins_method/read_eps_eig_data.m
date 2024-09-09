function data_out = read_eps_eig_data(run_in, label_in)
  % data_out = read_eigen_data(run_in, label_in)
  %
  % Reads the chart data from the previous solution [label_in] of [run_in],
  % grabs the unstable and stable eigenvectors and eigenvalues, and outputs
  % them in arrays.
  %
  % Input
  % ----------
  % run_in : str
  %     The string identifier for the previous COCO run that we will
  %     read information from.
  % label_in : int
  %     The integer solution label from the previous run.
  %
  % Output
  % ----------
  % data_out: data structure
  %     Data structure containing the previous run's values of eps1, esp2,
  %     theta, and the unstable and stable eigenvectors and eigenvalues.

  %----------------------------------%
  %     Read Data: Epsilon Stuff     %
  %----------------------------------%
  % Extract stored deviations of stable and unstable manifolds from
  % stationary equilibrium point
  [data, chart] = coco_read_solution('bcs_initial', run_in, label_in);
  epsilon = chart.x(data.epsilon_idx);

  % Normal vector
  normal = data.normal;
  % Intersection point
  pt0    = data.pt0;

  %--------------------------------%
  %     Read Data: Eigen Stuff     %
  %--------------------------------%
  % Unstable eigenvector and eigenvalue
  [data, chart] = coco_read_solution('bcs_eig_unstable', run_in, label_in);
  vu = chart.x(data.vu_idx);
  lu = chart.x(data.lu_idx);

  % Stable eigenvector and eigenvalue 1
  [data, chart] = coco_read_solution('bcs_eig_stable1', run_in, label_in);
  vs1 = chart.x(data.vs1_idx);
  ls1 = chart.x(data.ls1_idx);

  % Stable eigenvector and eigenvalue 2
  [data, chart] = coco_read_solution('bcs_eig_stable2', run_in, label_in);
  vs2 = chart.x(data.vs2_idx);
  ls2 = chart.x(data.ls2_idx);
  
  %----------------%
  %     Output     %
  %----------------%
  % Normal vector
  data_out.normal  = normal;
  % Intersection point
  data_out.pt0     = pt0;
  % Epsilon things
  data_out.epsilon = epsilon;

  % Eigenvectors
  data_out.vu      = vu;
  data_out.vs1     = vs1;
  data_out.vs2     = vs2;

  % Eigenvalues
  data_out.lu      = lu;
  data_out.ls1     = ls1;
  data_out.ls2     = ls2;

end