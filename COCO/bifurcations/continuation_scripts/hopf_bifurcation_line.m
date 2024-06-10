%-------------------------------------------------------------------------%
%%                         Hopf Bifurcation Loop                         %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.hopf_bifurcations;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'HB');
label_old = label_old(1);

% Print to console
fprintf('~~~ Third run (ode_HB2HB) ~~~ \n');
fprintf('Calculate line of Hopf bifurcation points H\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Set up COCO problem
prob = coco_prob();

% Set NAdapt to 1?
prob = coco_set(prob, 'cont', 'NAdapt', 5);

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2, 'h0', 1e-2, 'h_max', 1e-2);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set upper bound of continuation steps in each direction along solution
PtMX = 1000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);
% prob = coco_set(prob, 'cont', 'PtMX', [PtMX, 0]);
% prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set number of points
prob = coco_set(prob, 'coll', 'NTST', 700);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

% Continue from branching point
prob = ode_HB2HB(prob, '', run_old, label_old);

% Add saved point for double limit cycle calculation
prob = coco_add_event(prob, 'DL_PT', 'A', 6.715);

% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'}, p_range);
