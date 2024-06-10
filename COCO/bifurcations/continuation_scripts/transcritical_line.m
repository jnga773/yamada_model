%-------------------------------------------------------------------------%
%%                      Transcritical Points at A_T                      %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.transcritical;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = label_old(3);

% Print to console
fprintf('~~~ Fifth run (ode_ep2ep) ~~~ \n');
fprintf('Calculate line of transcritical points A_T\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Set up COCO problem
prob = coco_prob();

% Set NAdapt to 1?
prob = coco_set(prob, 'cont', 'NAdapt', 5);

% Set upper bound of continuation steps in each direction along solution
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

% Continue from branching point
prob = ode_ep2ep(prob, '', run_old, label_old);

% Run COCO continuation
coco(prob, run_new, [], 1, {'gamma', 'A'}, gamma_range);
