function save_stable_W_PO_data(run_in, label_in)
  % save_stable_W_PO_data(run_in, label_in)
  %
  % Reads periodic orbit solution data from COCO solution, calculates the
  % one-dimensional stable manifold of the "central" saddle point 'q', and
  % saves the data to './data/initial_PO.mat'.

  % Data matrix filename
  filename_out = './data_mat/stable_manifold_PO.mat';

  %-------------------%
  %     Read Data     %
  %-------------------%
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
  % Read periodic orbit solution
  data_out = load('./data_mat/initial_PO.mat');

  data_out.W1 = {W1_iso1, W1_iso2, W1_iso3};
  data_out.W2 = {W2_iso1, W2_iso2, W2_iso3};

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_out, '-struct', 'data_out');

end