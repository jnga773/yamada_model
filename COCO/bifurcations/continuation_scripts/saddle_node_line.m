%-------------------------------------------------------------------------%
%%                       Saddle-Node Points at A_S                       %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.saddle_nodes;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'SN');
label_old = label_old(1);

% Print to console
fprintf('~~~ Fourth run (ode_SN2SN) ~~~ \n');
fprintf('Calculate line of saddle-node points A_S\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Set up COCO problem
prob = coco_prob();

% Set NAdapt to 1?
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set upper bound of continuation steps in each direction along solution
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

% Continue from saddle-node
prob = ode_SN2SN(prob, '', run_old, label_old);

% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'}, p_range);
