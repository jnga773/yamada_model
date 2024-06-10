%-------------------------------------------------------------------------%
%%           Solve for Unstable Manifold towards data_lins.pt0            %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.lins_method.unstable_manifold;
% Which run this continuation continues from
run_old = run_names.approx_homo.continue_homoclinics;

% Grab the label for the previous run solution
label_old = 1;

% Print to console
fprintf("~~~ Lin's Method: First Run (ode_isol2coll) ~~~ \n");
fprintf('Continue unstable trajectory segment until we hit Sigma plane \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s\n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Initial Lin's Method data structure
data_lins = calc_lins_initial_conditions(run_old, label_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Construct instance of huxley continuation problem from initial data.
prob = coco_prob();

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set NTST size
prob = coco_set(prob, 'coll', 'NTST', 25);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set Continuation steps
PtMX = 400;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Construct first instance of 'coll' toolbox for unstable manifold
prob = ode_isol2coll(prob, 'unstable', funcs{:}, ...
                     data_lins.t0, data_lins.x_init_u, data_lins.pnames, data_lins.p0);
% Construct second instance of 'coll' toolbox for stable manifold
prob = ode_isol2coll(prob, 'stable', funcs{:}, ...
                     data_lins.t0, data_lins.x_init_s, data_lins.p0);

% Construct instance of 'ep' tool box to follow stationary point x_neg
% for initial conditions.
prob = ode_isol2ep(prob, 'x_neg', funcs{:}, data_lins.x_neg, ...
                   data_lins.p0);
% Construct instance of 'ep' tool box to follow stationary point x_pos
% for Lin's segment conditions.
prob = ode_isol2ep(prob, 'x_pos', funcs{:}, data_lins.x_pos, ...
                   data_lins.p0);

% Glue that shit together, haumi ;)
prob = glue_conditions(prob, data_lins);

% Run COCO
coco(prob, run_new, [], 1, {'seg_u', 'T1', 'T2'});
% coco(prob, run_new, [], 1, {'seg_u'});

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
% Grab maximum point of Sig_u
label_plot = coco_bd_labs(coco_bd_read(run_new), 'DelU');
label_plot = sort(label_plot);
label_plot = label_plot(2);

%--------------%
%     Plot     %
%--------------%
plot_homoclinic_manifold_run(run_new, label_plot, data_lins, 14, save_figure);

% plot_temporal_solution_single(run_new, label_plot, 15, save_figure);

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
