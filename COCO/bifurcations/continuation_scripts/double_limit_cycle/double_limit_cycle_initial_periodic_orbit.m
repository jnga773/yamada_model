%-------------------------------------------------------------------------%
%%               Solving for periodic solutions using ode45              %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.limit_cycle.initial_PO;
% Which run this continuation continues from
run_old = run_names.hopf_bifurcations;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'DL_PT');
label_old = max(label_old);

% Print to console
fprintf('~~~ Double Limit Cycle: Second Run (ode_SN2SN) ~~~ \n');
fprintf('Continue saddle-node of periodic orbits\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%---------------------------------------------%
%     Initial Periodic Orbit Continuation     %
%---------------------------------------------%
% Set up the COCO problem
prob = coco_prob();

% The value of 10 for 'NAdapt' implied that the trajectory discretisation
% is changed adaptively ten times before the solution is accepted.
prob = coco_set(prob, 'coll', 'NTST', 15);
prob = coco_set(prob, 'cont', 'NAdapt', 10);
% prob = coco_set(prob, 'po', 'bifus', 'off');
prob = coco_set(prob, 'coll', 'MXCL', false);

% Continue periodic orbit from initial solution
prob = ode_HB2po(prob, '', run_old, label_old);

% Fix the period
% prob = coco_xchg_pars(prob, 'gamma', 'po.period');

% Set step sizes
% h_temp = 1e-1;
% prob = coco_set(prob, 'cont', 'h_min', h_temp, 'h0', h_temp, 'h_max', h_temp);

% Set upper bound of continuation steps in each direction along solution
PtMX = 150;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

coco(prob, run_new, [], 1, {'gamma', 'A'}, gamma_range);

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
% plot_initial_periodic_orbits(t_PO, x_PO, run_new, save_figure);
% plot_initial_periodic_orbits_2D(t_PO, x_PO, run_new, save_figure);
% plot_initial_periodic_orbits_3D(t_PO, x_PO, run_new, save_figure);
