%-------------------------------------------------------------------------%
%%            Solve for Stable Manifold towards data_lins.pt0             %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.lins_method.stable_manifold;
% Which run this continuation continues from
run_old = run_names.lins_method.unstable_manifold;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'DelU');
label_old = sort(label_old);
label_old = label_old(2);

% Print to console
fprintf("~~~ Lin's Method: Second Run (ode_coll2coll) ~~~\n");
fprintf('Continue stable trajectory segment until we hit Sigma plane\n');
fprintf('Run name: %s\n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read epsilon and eig data from previous run
data_eps_eig = read_eps_eig_data(run_old, label_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Construct instance of huxley continuation problem from initial data.
prob = coco_prob();

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set Continuation steps
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue form previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'x_neg', run_old, label_old);
% Continue form previous 'ep' solution for the x_pos equilibrium point.
prob = ode_ep2ep(prob, 'x_pos', run_old, label_old);

% Glue that shit together, haumi ;)
prob = glue_conditions(prob, data_eps_eig);

% Run COCO
coco(prob, run_new, [], 1, {'seg_s', 'T1', 'T2'});
% coco(prob, run_new, [], 1, {'seg_s'});

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
% Grab maximum point of Sig_s
label_plot = coco_bd_labs(coco_bd_read(run_new), 'DelS');
label_plot = sort(label_plot);
label_plot = label_plot(1);

%--------------%
%     Plot     %
%--------------%
plot_homoclinic_manifold_run(run_new, label_plot, data_lins, 15, save_figure);

plot_temporal_solution_single(run_new, label_plot, 15, save_figure);

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
