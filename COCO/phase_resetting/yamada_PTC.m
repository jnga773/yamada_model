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

% % Hardcoded functions
% addpath('./functions/hardcoded/');
% addpath('./boundary_conditions/hardcoded/');
% Symbolic functions
addpath('./functions/symbolic/');
addpath('./boundary_conditions/symbolic/');

% Add continuation scripts
addpath('./continuation_scripts/');

% Add plotting scripts
addpath('./plotting_scripts/');

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

% Boundary conditions: Segments 1 and 2
% bcs_funcs.bcs_seg1_seg2 = {@bcs_PR_seg1_seg2};
bcs_funcs.bcs_seg1_seg2 = bcs_PR_seg1_seg2_symbolic();

% Boundary conditions: Segment 3
% bcs_funcs.bcs_seg3 = {@bcs_PR_seg3};
bcs_funcs.bcs_seg3 = bcs_PR_seg3_symbolic();

% Boundary conditions: Segment 4
% bcs_funcs.bcs_seg4 = {@bcs_PR_seg4};
bcs_funcs.bcs_seg4 = bcs_PR_seg4_symbolic();

%-------------------------------------------------------------------------%
%%                   Increasing Pertubation Amplitude                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.phase_reset_perturbation = 'run09_phase_reset_perturbation';
% Which run this continuation continues from
run_old = 'run08_compute_floquet_bundle_2_w';

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'NORM1');
label_old = label_old(1);

% Print to console
fprintf("~~~ Phase Transition Curve: First Run ~~~ \n");
fprintf('Increase perturbation amplitude \n');
fprintf('Run name: %s \n', run_names.phase_reset_perturbation);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read data from _previous solution
data_PR = calc_PR_initial_conditions(run_old, label_old);

% Set some parameters
% k: Integer number of periods
data_PR.p0(data_PR.p_maps.k) = 25;
% theta_perturb: Angle at which perturbation is applied in the (G-I) plane
data_PR.p0(data_PR.p_maps.theta_perturb) = 0.0;
% data_PR.p0(data_PR.p_maps.theta_perturb) = 0.5 * pi;

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set tolerance
prob = coco_set(prob, 'corr', 'TOL', 5e-7);

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 5e-5);
prob = coco_set(prob, 'cont', 'h0', 1e-3);
prob = coco_set(prob, 'cont', 'h_max', 1e2);

% Set adaptive mesh
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set number of steps
PtMX = 1000;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

%------------------%
%     Set NTST     %
%------------------%
% NTST values
NTST(1) = 25;
NTST(2) = 20;
NTST(3) = 20;
NTST(4) = 25 * data_PR.p0(data_PR.p_maps.k);

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST(1));
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST(2));
prob = coco_set(prob, 'seg3.coll', 'NTST', NTST(3));
prob = coco_set(prob, 'seg4.coll', 'NTST', NTST(4));

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

% Set norm to int
prob = coco_set(prob, 'norm', 'Inf');

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply all boundary conditions, glue parameters together, and
% all that other good COCO stuff. Looking the function file
% if you need to know more ;)
prob = apply_PR_boundary_conditions(prob, data_PR, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Array of values for special event
SP_values = [linspace(0.0, 0.25, 21), linspace(0.30, 2.0, 41)];

% When the parameter we want (from param) equals a value in A_vec
prob = coco_add_event(prob, 'SP', 'A_perturb', SP_values);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[-1e-4, max(SP_values)], [], [], [0.99, 1.01], []};
coco(prob, run_names.phase_reset_perturbation, [], 1, {'A_perturb', 'theta_new', 'eta', 'mu_s', 'T'}, prange);

%--------------------%
%     Test Plots     %
%--------------------%
% Label of solution to plot
label_plot = sort(coco_bd_labs(coco_bd_read(run_names.phase_reset_perturbation), 'SP'));
label_plot = label_plot(end);

% Plot some stuff my g
plot_phase_reset_phase_space(run_names.phase_reset_perturbation, label_plot, 1);

% Plot perturbation amplitude against theta_new
plot_A_perturb_theta_new(run_names.phase_reset_perturbation);

%-------------------------------------------------------------------------%
%%                   Phase Response Curve Calculation                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.phase_transition_curve = 'run10_phase_reset_PTC_single';
% Which run this continuation continues from
run_old = run_names.phase_reset_perturbation;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'SP');
label_old = label_old(end);

% Print to console
fprintf("~~~ Phase Reset: Second Run (phase_reset_2_theta_perturb.m) ~~~ \n");
fprintf('Fix A_perturb and continue in theta_perturb \n');
fprintf('Run name: %s \n', run_names.phase_transition_curve);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set tolerance
prob = coco_set(prob, 'corr', 'TOL', 5e-7);

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 5e-2);
prob = coco_set(prob, 'cont', 'h0', 1e-1);
prob = coco_set(prob, 'cont', 'h_max', 1e2);

% Set adaptive mesh
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set number of steps
prob = coco_set(prob, 'cont', 'PtMX', 750);

%-------------------------------------------%
%     Continue from Trajectory Segments     %
%-------------------------------------------%
% Segment 1
prob = ode_coll2coll(prob, 'seg1', run_old, label_old);
% Segment 2
prob = ode_coll2coll(prob, 'seg2', run_old, label_old);
% Segment 3
prob = ode_coll2coll(prob, 'seg3', run_old, label_old);
% Segment 4
prob = ode_coll2coll(prob, 'seg4', run_old, label_old); 

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply all boundary conditions, glue parameters together, and
% all that other good COCO stuff. Looking the function file
% if you need to know more ;)
prob = apply_PR_boundary_conditions(prob, data_PR, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Run COCO continuation
prange = {[0.0, 1.0], [], [], [0.99, 1.01], [], []};
coco(prob, run_names.phase_transition_curve, [], 1, {'theta_old', 'theta_new', 'eta', 'mu_s', 'T', 'A_perturb'}, prange);

%--------------------%
%     Test Plots     %
%--------------------%
% Plot first SP solution
plot_phase_reset_phase_space(run_names.phase_transition_curve, label_plot(1), 1);

% Plot phase transition curve (PTC)
plot_phase_transition_curve(run_names.phase_transition_curve);

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%