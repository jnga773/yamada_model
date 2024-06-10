function data_out = read_lingap_data(run_in, label_in)
  % data_out = read_lingap_data(run_in, label_in)
  %
  % Reads the chart data from the previous solution [label_in] of [run_in],
  % grabs the Lin gap and vectors data.
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
  %     Data structure containing the previous run's values for the Lin gap
  %     vector and distance.

  %----------------------------------%
  %     Read Data: Epsilon Stuff     %
  %----------------------------------%
  % Read Linsgap data structure
  data_lingap = coco_read_solution('lins_data', run_in, label_in);

  % Extract stored lingap value from previous run
  [data, chart] = coco_read_solution('bcs_lingap', run_in, label_in);
  lingap = chart.x(data.lingap_idx);
  
  %----------------%
  %     Output     %
  %----------------%
  % Data structure
  data_out        = data_lingap;
  % Lingap value
  data_out.lingap = lingap;

end