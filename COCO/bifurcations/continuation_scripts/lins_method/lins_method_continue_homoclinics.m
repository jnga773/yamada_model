%-------------------------------------------------------------------------%
%%                     Parametrise the Heteroclinic                      %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.lins_method.continue_homoclinics;
% Which run this continuation continues from
% run_old = run_names.lins_method.close_lingap;
run_old = run_names.lins_method.close_eps2;

% Label for previous run solution
% label_old = coco_bd_labs(coco_bd_read(run_old), 'Lin0');
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = max(label_old);

% Print to console
fprintf("~~~ Lin's Method: Sixth Run (ode_coll2coll) ~~~ \n");
fprintf('Continue constrained segments to find parametrisation of homoclinic \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read epsilon and eig data from previous run
data_eps_eig = read_eps_eig_data(run_old, label_old);

% Calculate Lin gap and vector
data_lingap = read_lingap_data(run_old, label_old);

%-------------------------------------------%
%     Initialise Problem Structure (OLD)    %
%-------------------------------------------%
% Construct instance of huxley continuation problem from initial data.
prob = coco_prob();

% Set NTST size
% prob = coco_set(prob, 'coll', 'NTST', 100);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set step sizes
% h = 1e-3;
% prob = coco_set(prob, 'cont', 'h_min', h, 'h0', h, 'h_max', h);

% Set frequency of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 5);

% Set Continuation steps
PtMX = 500;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);
% prob = coco_set(prob, 'cont', 'PtMX', [PtMX, 0]);
% prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

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

% Apply Lin's conditions
prob = glue_lingap_conditions(prob, data_lingap);

% Calculate Bogdanov Takens Point
% gamma_BT = -(1 + (B * (a - 1)) - 2 * sqrt(B * (a - 1))) / (sqrt(B * (a - 1)) * (1 - a - sqrt(B * (a - 1))));

% Run COCO
p_range_test = {A_range, ...
                [0.0, -(1 + (B * (a - 1)) - 2 * sqrt(B * (a - 1))) / (sqrt(B * (a - 1)) * (1 - a - sqrt(B * (a - 1))))]};
coco(prob, run_new, [], 1, {'A', 'gamma', 'eps1', 'eps2', 'theta', 'seg_u'}, p_range_test);
% coco(prob, run_new, [], 1, {'A', 'gamma', 'T2', 'eps1', 'theta', 'seg_u'}, p_range);

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
% Choose solution label
% label_plot = coco_bd_labs(coco_bd_read(run_new), 'EP');
% label_plot = sort(label_plot);
% label_plot = label_plot(2);
label_plot = 6;

% Read solution
[solu, ~] = coll_read_solution('unstable', run_new, label_plot);
[sols, ~] = coll_read_solution('stable', run_new, label_plot);
xu_plot = solu.xbp;
xs_plot = sols.xbp;

%--------------%
%     Plot     %
%--------------%
plot_homoclinic_manifold_run(run_new, label_plot, data_lins, 19, save_figure);

plot_temporal_solution_single(run_new, label_plot, 20, save_figure);

plot_temporal_solutions(run_new, 21, save_figure);

compare_homoclinic_bifurcations(run_names, save_figure);

%--------------------------%
%     Print to Console     %
%--------------------------%
[sol1, ~] = coll_read_solution('unstable', run_new, label_plot);
[sol2, ~] = coll_read_solution('stable', run_new, label_plot);
[sol_neg, ~] = ep_read_solution('x_neg', run_new, label_plot);

fprintf('Print Start and End Points to Console\n');
fprintf('Equilibrium point       = (%.3f, %.3f, %.3f)\n', sol_neg.x);
fprintf('Unstable starting point = (%.3f, %.3f, %.3f)\n', sol1.xbp(1, :));
fprintf('Unstable ending point   = (%.3f, %.3f, %.3f)\n', sol1.xbp(end, :));
fprintf('Stable starting point   = (%.3f, %.3f, %.3f)\n', sol2.xbp(1, :));
fprintf('Stable ending point     = (%.3f, %.3f, %.3f)\n', sol2.xbp(end, :));

