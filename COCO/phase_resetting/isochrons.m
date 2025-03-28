%=========================================================================%
%                     YAMADA MODEL (Phase Resetting)                      %
%=========================================================================%
% We compute the phase resetting of an attracting periodic orbit of the
% Yamada model:
%                     G' = \gamma (A - G - G I) ,
%                     Q' = \gamma (B - Q - a Q I) ,
%                     I' = (G - Q - 1) I ,
% where G is the gain, Q is the absorption, and I is the intensity of the
% laser. The system is dependent on four parameters: the pump current on
% the gain, A (or A); the relative absoprtion, B and a; and the decay
% time of the gain, \gamma.

% Clear plots
close('all');

% Clear workspace
clear;
clc;

% Add equation/functions to path
addpath('./functions/');
% Add field functions to path
addpath('./functions/fields/');
% Add boundary condition functions to path
addpath('./functions/bcs/');
% Add SymCOCO files to path
addpath('./functions/symcoco/');

% Add continuation scripts
addpath('./continuation_scripts/phase_reset/');
% Add plotting scripts
addpath('./plotting_scripts/isochrons/');

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Phase Reset Segment 1: Functions
% func.seg1 = {@func_seg1};
funcs.seg1 = func_seg1_symbolic();

% Phase Reset: Segment 2
% funcs.seg2 = {@func_seg2};
funcs.seg2 = func_seg2_symbolic();

% Phase Reset: Segment 3
% funcs.seg3 = {@func_seg3};
funcs.seg3 = func_seg3_symbolic();

% Phase Reset: Segment 4
% funcs.seg4 = {@func_seg4};
funcs.seg4 = func_seg4_symbolic();

% Boundary conditions: Period
% bcs_funcs.bcs_T = {@bcs_T};
bcs_funcs.bcs_T = bcs_T_symbolic();

% Boundary conditions: Phase-resetting segments
% bcs_funcs.bcs_PR = {@bcs_isochron};
bcs_funcs.bcs_PR = bcs_isochron_symbolic();

%-------------------------------------------------------------------------%
%%            Move Along Periodic Orbit (theta_old, theta_new)           %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.isochron_initial = 'run01_isochron_initial';
run_new = run_names.isochron_initial;

% Print to console
fprintf("~~~ Isochrons: First Run ~~~ \n");
fprintf('Continue around PO with zero perturbation \n');
fprintf('Run name: %s \n', run_new);

%-------------------%
%     Read Data     %
%-------------------%
% Set periodicity
k = 20;

% Set initial conditions from previous solutions
data_PR = calc_initial_solution_PR('./data_mat/solution_VAR.mat', k, isochron=true);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
% prob = coco_set(prob, 'cont', 'h_min', 5e-1);
% prob = coco_set(prob, 'cont', 'h0', 1e0);
prob = coco_set(prob, 'cont', 'h_max', 1e1);

% Set adaptive mesh
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set number of steps
PtMX = 1000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

% Set norm to int
prob = coco_set(prob, 'cont', 'norm', inf);

% Set MaxRes and al_max
% prob = coco_set(prob, 'cont', 'MaxRes', 10);
% prob = coco_set(prob, 'cont', 'al_max', 25);

%------------------%
%     Set NTST     %
%------------------%
% In calc_PR_initial conditions, we define segment 4 as having 'k' periods,
% where 'k' is an integer. This is the perturbed segment, that may have to
% orbit the unperturbed periodic orbit many times before "resetting". Hence
% we have set the NTST for this segment (NTST(4)) as k * 50.
NTST(1) = 25;
NTST(2) = 10;
NTST(3) = 10;
NTST(4) = 25 * k;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST(1));
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST(2));
prob = coco_set(prob, 'seg3.coll', 'NTST', NTST(3));
prob = coco_set(prob, 'seg4.coll', 'NTST', NTST(4));

% % Set min and max NTST for segment 4
% prob = coco_set(prob, 'seg4.coll', 'NTSTMN', 10*k);
% prob = coco_set(prob, 'seg4.coll', 'NTSTMX', 100*k);

%------------------------------------%
%     Create Trajectory Segments     %
%------------------------------------%
% Here we set up the four segments to calculate the phase resetting curve:
% Segment 1 - Trajectory segment of the periodic orbit from the zero-phase
%             point (gamma_0) to the point where the perturbed trajectory 
%             comes close to the periodic orbit (at theta_new).
prob = ode_isol2coll(prob, 'seg1', funcs.seg1{:}, ...
                     data_PR.t_seg1, data_PR.x_seg1, data_PR.p0);

% Segment 2 - Trajectory segment of the periodic from the end of Segment 1
%             (at theta_new) back to the zero-phase point (gamma_0).
prob = ode_isol2coll(prob, 'seg2', funcs.seg2{:}, ...
                     data_PR.t_seg2, data_PR.x_seg2, data_PR.p0);

% Segment 3 - Trajectory segment of the periodic orbit from the zero-phase
%             point (gamma_0) to the point at which the perturbation is
%             applied (theta_old).
prob = ode_isol2coll(prob, 'seg3', funcs.seg3{:}, ...
                     data_PR.t_seg3, data_PR.x_seg3, data_PR.p0);   

% Segment 4 - Trajectory segment of the perturbed trajectory, from
%             theta_old to theta_new.
prob = ode_isol2coll(prob, 'seg4', funcs.seg4{:}, ...
                     data_PR.t_seg4, data_PR.x_seg4, data_PR.p0);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply all boundary conditions, glue parameters together, and
% all that other good COCO stuff. Looking the function file
% if you need to know more ;)
prob = apply_boundary_conditions_PR(prob, data_PR, bcs_funcs, isochron=true);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Array of values for special event
% SP_values = [linspace(0.0, 0.25, 21), linspace(0.30, 2.0, 21)];
% SP_values = 0.0 : 0.05 : 1.0;

% Save solutions at zero-phase point, where I is max, and where I is min
SP_values = [0.0, 0.1858, 0.6768];

% When the parameter we want (from param) equals a value in A_vec
prob = coco_add_event(prob, 'SP', 'theta_old', SP_values);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
% prange = {[0.0, 1.0], [0.0, 1.0], [], [0.99, 1.01], []};
prange = {[0.0, 1.0], [0.0, 1.0], [], [0.99, 1.01], []};
bdtest = coco(prob, run_new, [], 1, ...
              {'theta_old', 'theta_new', 'eta', 'mu_s', 'T'}, prange);

%--------------------%
%     Test Plots     %
%--------------------%
% Label of solution to plot
label_plot = sort(coco_bd_labs(coco_bd_read(run_new), 'SP'));
label_plot = label_plot(2);

% Plot first SP solution
plot_phase_reset_phase_space(run_new, label_plot(1), 1);

%-------------------------------------------------------------------------%
%%                  Calculate Isochrons in d_y Plane                     %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.single_isochron_dy = 'run02_single_isochron_test';
run_new = run_names.single_isochron_dy;
% Which run this continuation continues from
run_old = run_names.isochron_initial;

% Continuation point
% label_old = sort(coco_bd_labs(coco_bd_read(run_old), 'SP'));
% label_old = label_old(1);
label_old = 1;

% Print to console
fprintf("~~~ Isochrons: Second Run ~~~ \n");
fprintf('Fix theta_old and theta_new, and follow isochrons in one direciton \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
run_isochron_continuation(this_run_name, run_old, label_old, data_PR, bcs_funcs, {'d_x', 'd_y'});

%--------------------%
%     Test Plots     %
%--------------------%
% Plot phase transition curve (PTC)
plot_single_isochron(run_new);

%-------------------------------------------------------------------------%
%%                Calculate Isochrons in (d_x, d_z) Plane                %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.single_isochron = 'run03_single_isochron_test';
run_new = run_names.single_isochron;
% Which run this continuation continues from
run_old = run_names.single_isochron_dy;

% Continuation point
label_old = sort(coco_bd_labs(coco_bd_read(run_old), 'SP'));
label_old = label_old(3);

% Print to console
fprintf("~~~ Isochrons: Third Run ~~~ \n");
fprintf('Continue isochron in (dx, dz) plane \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
run_isochron_continuation(run_new, run_old, label_old, data_PR, bcs_funcs, {'d_x', 'd_z'});
prange = {[], [], [], [0.99, 1.01], [], [-4, 6], [-4, 6], []};
bdtest = coco(prob, run_new, [], 1, {'d_x', 'd_z', 'eta', 'mu_s', 'T', 'iso1', 'iso2', 'iso3'}, prange);

%--------------------%
%     Test Plots     %
%--------------------%
% Plot phase transition curve (PTC)
plot_single_isochron(run_new);

%-------------------------------------------------------------------------%
%%                       Isochron Scan - Multiple                        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.isochron_scan = 'run04_isochron_scan';
run_new = run_names.isochron_scan;
% Which run this continuation continues from
run_old = run_names.single_isochron_dy;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'SP');

% Print to console
fprintf("~~~ Isochron: Fourth Run ~~~ \n");
fprintf('Calculate slices of the isochrons \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from SP points in run: %s \n', run_old);

%---------------------------------%
%     Cycle through SP labels     %
%---------------------------------%
% Set number of threads
M = 7;
parfor (run = 1 : length(label_old), M)
  % Label for this run
  this_run_label = label_old(run);

  % Data directory for this run
  fprintf('\n Continuing from point %d in run: %s \n', this_run_label, run_old);

  this_run_name = {run_new; sprintf('run_%02d', run)};

  % Run continuation
  isochron_scan(this_run_name, run_old, this_run_label, data_PR, bcs_funcs, {'d_x', 'd_y'});

end

%--------------------%
%     Test Plots     %
%--------------------%
% Plot isochron_scan
plot_phase_reset_phase_space({run_new, 'run_08'}, 5, 1);

plot_single_isochron({run_new, 'run_08'});
plot_isochron_scan(run_new);

% Save isochron save data
save_isochron_data(run_new, './data_mat/isochron_scan_data_COCO.mat');

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%