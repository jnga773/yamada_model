%-------------------------------------------------------------------------%
%%                   Lin's Method: Closing the Lin Gap                   %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.lins_method.close_lingap;
% Which run this continuation continues from
run_old = run_names.lins_method.stable_manifold;
% run_old = run_names.lins_method.close_eps2;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'DelS');
label_old = sort(label_old);
label_old = label_old(1);
% label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
% label_old = max(label_old);

% Print to console
fprintf("~~~ Lin's Method: Third Run (ode_coll2coll) ~~~ \n");
fprintf('Close the Lin Gap on the Sigma Plane \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read epsilon and eig data from previous run
data_eps_eig = read_eps_eig_data(run_old, label_old);

% Calculate Lin gap and vector
data_lingap = calculate_lingap_vector(run_old, label_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Construct instance of huxley continuation problem from initial data.
prob = coco_prob();

% Set NAdapt
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set Continuation steps
PtMX = 1500;
% prob = coco_set(prob, 'cont', 'PtMX', PtMX);
prob = coco_set(prob, 'cont', 'PtMX', [PtMX, 0]);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue form previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'x_neg', run_old, label_old);
% Continue form previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'x_pos', run_old, label_old);

% Glue that shit together, haumi ;)
prob = glue_conditions(prob, data_eps_eig);

% Initial lingap
lingap = data_lingap.lingap0;

% Apply Lin's conditions
prob = glue_lingap_conditions(prob, data_lingap);

% Assign 'gamma' to the set of active continuation parameters and 'po.period'
% to the set of inactive continuation parameters, thus ensuring that the
% latter is fixed during root finding.
% prob = coco_xchg_pars(prob, 'T1', 'A');

% Run COCO
% coco(prob, run_new, [], 1, {'lingap', 'T1', 'T2', 'theta', 'A', 'seg_s'}, [0, 1.5*lingap]);
coco(prob, run_new, [], 1, {'lingap', 'eps1', 'eps2', 'theta', 'gamma', 'seg_u'}, [0, 1.5*lingap]);

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
% Find good label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'Lin0');
label_plot = sort(label_plot);
label_plot = label_plot(1);

%--------------%
%     Plot     %
%--------------%
plot_homoclinic_manifold_run(run_new, label_plot, data_lins, 16, save_figure);

plot_temporal_solution_single(run_new, label_plot, 17, save_figure);

%--------------------------%
%     Print to Console     %
%--------------------------%
[sol1, ~] = coll_read_solution('unstable', run_new, label_plot);
[sol2, ~] = coll_read_solution('stable', run_new, label_plot);
[sol_neg, ~] = ep_read_solution('x_neg', run_new, label_plot);
[sol_pos, ~] = ep_read_solution('x_pos', run_new, label_plot);

fprintf('Print Start and End Points to Console\n');
fprintf('Equilibrium point (neg) = (%.3f, %.3f, %.3f)\n', sol_neg.x);
fprintf('Equilibrium point (pos) = (%.3f, %.3f, %.3f)\n', sol_pos.x);
fprintf('Unstable starting point = (%.3f, %.3f, %.3f)\n', sol1.xbp(1, :));
fprintf('Unstable ending point   = (%.3f, %.3f, %.3f)\n', sol1.xbp(end, :));
fprintf('Stable starting point   = (%.3f, %.3f, %.3f)\n', sol2.xbp(1, :));
fprintf('Stable ending point     = (%.3f, %.3f, %.3f)\n', sol2.xbp(end, :));

%-------------------------------------------------------------------------%
%%                               FUNCTIONS                               %%
%-------------------------------------------------------------------------%
function data_out = calculate_lingap_vector(run_in, label_in)
  % data_out = calculate_lingap_vector(run_in)
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
  
  % Add normalised Lin gap vector to data_out
  data_lins.vgap = vgap / vgap_norm;

  % Add value of norm to data_out
  data_lins.lingap0 = vgap_norm;
  data_lins.lingap  = vgap_norm;

  %----------------%
  %     Output     %
  %----------------%
  data_out      = data_lins;
  data_out.xdim = datau.xdim;
  data_out.pdim = datau.pdim;

end
