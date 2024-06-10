%-------------------------------------------------------------------------%
%%                           Move Hopf A Value                           %%
%-------------------------------------------------------------------------%
% Continuing from a Hopf bifurcation with 'ode_HB2HB', we vary
% the 'A' parameter to A = 6.4

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.move_hopf;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'HB');
label_old = label_old(1);

% Print to console
fprintf("~~~ Initial Periodic Orbit: Third Run (move_hopf.m) ~~~ \n");
fprintf('Move the gamma value \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2, 'h0', 1e-2, 'h_max', 1e-2);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 25);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

% Set upper bound of continuation steps in each direction along solution
PtMX = 300;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Initial solution to periodic orbit (COLL Toolbox)
prob = ode_HB2HB(prob, '', run_old, label_old);

% Even for Nonlinear Photonics abstract
prob = coco_add_event(prob, 'H_PT', 'A', 7.3757);

% Set event for A = 7.5
% prob = coco_add_event(prob, 'H_PT', 'A', 6.9);

% Run COCO
coco(prob, run_new, [], 1, {'A', 'gamma'}, A_range);
