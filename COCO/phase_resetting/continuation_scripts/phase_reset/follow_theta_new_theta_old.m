%-------------------------------------------------------------------------%
%%                   Phase Response Curve Calculation                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.phase_reset_thetas;
% Which run this continuation continues from
run_old = run_names.compute_floquet_2;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'NORM1');
label_old = label_old(1);

% Print to console
fprintf("~~~ Phase Reset: First Run (phase_reset_1_perturbation.m) ~~~ \n");
fprintf('Calculate phase resetting curve \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read data from _previous solution
data_PR = calc_PR_initial_conditions(run_old, label_old);

% Set some parameters

% k: Integer number of periods
data_PR.p0(data_PR.p_maps.k) = 25;
% theta_perturb: Angle at which perturbation is applied
data_PR.p0(data_PR.p_maps.theta_perturb) = 0.0;
% phi_perturb: Azimuthal angle at which perturbation is applied
% data_PR.p0(data_PR.p_maps.phi_perturb) = 0.0;
data_PR.p0(data_PR.p_maps.phi_perturb) = 0.5 * pi;

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
PtMX = 800;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% NTST values
NTST(1) = 50;
NTST(2) = 20;
NTST(3) = 20;
NTST(4) = 50 * data_PR.p0(data_PR.p_maps.k);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

%------------------------------------%
%     Create Trajectory Segments     %
%------------------------------------%
% Here we set up the four segments to calculate the phase resetting curve:
% Segment 1 - Trajectory segment of the periodic orbit from the zero-phase
%             point (gamma_0) to the point where the perturbed trajectory 
%             comes close to the periodic orbit (at theta_new).
prob = coco_set(prob, 'coll', 'NTST', NTST(1));
prob = ode_isol2coll(prob, 'seg1', funcs.seg1{:}, ...
                     data_PR.t_seg1, data_PR.x_seg1, data_PR.p0);

% Segment 2 - Trajectory segment of the periodic from the end of Segment 1
%             (at theta_new) back to the zero-phase point (gamma_0).
prob = coco_set(prob, 'coll', 'NTST', NTST(2));
prob = ode_isol2coll(prob, 'seg2', funcs.seg2{:}, ...
                     data_PR.t_seg2, data_PR.x_seg2, data_PR.p0);

% Segment 3 - Trajectory segment of the periodic orbit from the zero-phase
%             point (gamma_0) to the point at which the perturbation is
%             applied (theta_old).
prob = coco_set(prob, 'coll', 'NTST', NTST(3));
prob = ode_isol2coll(prob, 'seg3', funcs.seg3{:}, ...
                     data_PR.t_seg3, data_PR.x_seg3, data_PR.p0);   

% Segment 4 - Trajectory segment of the perturbed trajectory, from
%             theta_old to theta_new.
prob = coco_set(prob, 'coll', 'NTST', NTST(4));
prob = ode_isol2coll(prob, 'seg4', funcs.seg4{:}, ...
                     data_PR.t_seg4, data_PR.x_seg4, data_PR.p0);       
                     
%----------------------------%
%     Create EP Instance     %
%----------------------------%
% We also want to follow the equilibrium point / singularity in the centre
% of the periodic orbit, so we do so with ode_isol2ep.
prob = ode_isol2ep(prob, 'singularity', funcs.field{:}, ...
                   data_PR.xpos, data_PR.p0(1:pdim));

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply all boundary conditions, glue parameters together, and
% all that other good COCO stuff. Looking the function file
% if you need to know more ;)
prob = glue_PR_conditions(prob, data_PR, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Array of values for special event
SP_values = 0.0 : 0.01 : 2.0;

% When the parameter we want (from param) equals a value in A_vec
prob = coco_add_event(prob, 'SP', 'theta_old', SP_values);

% Run COCO continuation
prange = {[0.0, 2.0], [0.0, 2.0], [-1e-4, 1e-2], [0.99, 1.01]};
coco(prob, run_new, [], 1, {'theta_old', 'theta_new', 'eta', 'mu_s'}, prange);

%-------------------------------------------------------------------------%
%%                            Testing Things                             %%
%-------------------------------------------------------------------------%
% Label of solution to plot
label_plot = sort(coco_bd_labs(coco_bd_read(run_new), 'SP'));
label_plot = label_plot(4);
% label_plot = 1;

% Plot some stuff my g
plot_phase_reset_phase_space(run_new, label_plot, 2, save_figure, true);

% Plot perturbation amplitude against theta_new
plot_A_perturb_theta_new(run_new, save_figure);

% Plot big pic
% plot_phase_reset_theta_new(run_new, save_figure);
