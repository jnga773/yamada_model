%-------------------------------------------------------------------------%
%%                 Finding Periodic Orbit Folding Point                  %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.limit_cycle.follow_limit_cycle;
% Which run this continuation continues from
run_old = run_names.limit_cycle.initial_PO;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'SN');
label_old = max(label_old);

% Print to console
fprintf('~~~ Double Limit Cycle: Second Run (ode_SN2SN) ~~~ \n');
fprintf('Continue saddle-node of periodic orbits\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Set up the COCO problem
prob = coco_prob();

prob = ode_SN2SN(prob, '', run_old, label_old);

% The value of 10 for 'NAdapt' implied that the trajectory discretisation
% is changed adaptively ten times before the solution is accepted.
% prob = coco_set(prob, 'coll', 'NTST', 15);
% prob = coco_set(prob, 'cont', 'NAdapt', 30);
% prob = coco_set(prob, 'po', 'bifus', 'off');
% prob = coco_set(prob, 'coll', 'MXCL', false);

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 2.5e0, 'h0', 2.5e0, 'h_max', 2.5e0);

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Set upper bound of continuation steps in each direction along solution
PtMX =  5000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set NPR to save every 100 steps
prob = coco_set(prob, 'cont', 'NPR', 100);

% Continue periodic orbit from initial solution
coco(prob, run_new, [], 1, {'A', 'gamma'}, p_range);
