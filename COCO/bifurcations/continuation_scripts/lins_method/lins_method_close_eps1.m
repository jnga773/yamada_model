%-------------------------------------------------------------------------%
%%                        Close the Distance eps1                        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.lins_method.close_eps1;
% Which run this continuation continues from
run_old = run_names.lins_method.close_lingap;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'Lin0');
label_old = sort(label_old);
label_old = label_old(1);

% Print to console
fprintf("~~~ Lin's Method: Fourth Run (ode_coll2coll) ~~~ \n");
fprintf('Close epsilon gap until eps1=1e-8 \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read epsilon and eig data from previous run
data_eps_eig = read_eps_eig_data(run_old, label_old);

% Calculate Lin gap and vector
data_lingap = read_lingap_data(run_old, label_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Construct instance of huxley continuation problem from initial data.
prob = coco_prob();

% Set NAdapt
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set step sizes
h = 5e-1;
prob = coco_set(prob, 'cont', 'h_min', h, 'h0', h, 'h_max', h);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set Continuation steps
PtMX = 400;
% prob = coco_set(prob, 'cont', 'PtMX', PtMX);
prob = coco_set(prob, 'cont', 'PtMX', [PtMX, 0]);

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

% Add event for when eps1 gets small enough
prob = coco_add_event(prob, 'EPS1', 'eps1', 1e-5);

% Run COCO
coco(prob, run_new, [], 1, {'eps1', 'theta', 'gamma', 'T1', 'T2', 'seg_u'});
% coco(prob, run_new, [], 1, {'eps1', 'theta', 'seg_u', 'gamma'});

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
% Find good label to plot
% label_plot = coco_bd_labs(coco_bd_read(run_new), 'EPS1');
% label_plot = label_plot(1);
label_plot = coco_bd_labs(coco_bd_read(run_new), 'EP');
label_plot = sort(label_plot);
label_plot = label_plot(end);

%--------------%
%     Plot     %
%--------------%
plot_homoclinic_manifold_run(run_new, label_plot, data_lins, 17, save_figure);

plot_temporal_solution_single(run_new, label_plot, 18, save_figure);

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
