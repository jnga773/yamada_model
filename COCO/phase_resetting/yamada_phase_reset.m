% Clear plots
close('all');

% Clear workspace
clear;
clc;

% Add equation/functions to path
addpath('./functions/');
addpath('./boundary_conditions/');
addpath('./continuation_scripts/');
addpath('./plotting_scripts/');

% Symbolic functions
addpath('./functions/symbolic/');
addpath('./boundary_conditions/symbolic/');

% Save figures switch
% save_figure = true;
save_figure = false;

%-------------------------------------------------------------------------%
%%                    YAMADA MODEL (Phase Resetting)                     %%
%-------------------------------------------------------------------------%
% We compute the phase resetting of an attracting periodic orbit of the
% Yamada model:
%     G' = \gamma (A - G - G I) ,
%     Q' = \gamma (B - Q - a Q I) ,
%     I' = (G - Q - 1) I ,
% where G is the gain, Q is the absorption, and I is the intensity of the
% laser. The system is dependent on four parameters: the pump current on
% the gain, A (or A); the relative absoprtion, B and a; and the decay
% time of the gain, \gamma.

% Here we compute isochrons of the limit cycles present in the Yamada
% model, and plot phase resetting curves.

%--------------------%
%     Parameters     %
%--------------------%
% Because we will only be looking at the (A, \gamma) plane, we will be
% setting values for a and B.
B = 5.8;
a = 1.8;

% Set some initial values for \gamma and A
gamma = 0.10;
A = 6.6;

%--------------------%
%     COCO Setup     %
%--------------------%
% Parameter names
pnames = {'gamma', 'A', 'B', 'a'};

% Initial parameter values
p0 = [gamma; A; B; a];

% Initial state values
x0 = [A; B; 0];

% Parameter ranges
gamma_range = [0.0, 0.25];
A_range = [5.0, 11.0];
p_range = {A_range, gamma_range};

% State dimensions
pdim = length(p0);
xdim = length(x0);

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Vector field: Functions
% funcs.field = {@yamada, @yamada_DFDX, @yamada_DFDP};
funcs.field = yamada_symbolic();

% Boundary conditions: Periodic orbit
% bcs_funcs.bcs_PO = {@bcs_PO};
bcs_funcs.bcs_PO = bcs_PO_symbolic();

% Boundary conditions: Period
% bcs_funcs.bcs_T = {@bcs_T};
bcs_funcs.bcs_T = bcs_T_symbolic();

%-------------------------------------------------------------------------%
%%                         Initial Continuation                          %%
%-------------------------------------------------------------------------%
% We set up the continuation problem by first continuing the equilibrium
% point x0. The equilibrium point undergoes a Hopf bifurcation, from which
% a family of periodic orbits originate.

% Add continuation scripts to path
addpath('./continuation_scripts/initial_periodic_orbit/');
addpath('./plotting_scripts/initial_periodic_orbit/');

%-------------------------------------%
%%     Compute Equilibrium Point     %%
%-------------------------------------%
% We compute and continue the equilibrium point of the model using
% the 'EP' toolbox constructor 'ode_isol2ep'.
run_names.initial_EP = 'run01_initial_EP';

% Run continuation script
initial_equilibrium_point;

%----------------------------------------%
%%     Continue from Branching Point     %
%----------------------------------------%
% Continue from the branching point to find a Hopf bifurcation
run_names.branching_point = 'run02_branching_point_continuation';

% Run continuation script
equilibrium_from_branching_point;

%-----------------------------------------%
%%     Move A and \gamma on Hopf Line     %
%-----------------------------------------%
% Continuing from a Hopf bifurcation with 'ode_HB2HB', we vary 'A' and 
% 'gamma' until we're kind of far from the Homoclinic.
run_names.move_hopf = 'run03_move_hopf';

% Run continuation script
move_hopf;

%----------------------------------%
%%     Hopf to Periodic Orbit     %%
%----------------------------------%
% Continue a family of periodic orbits emanating from the Hopf
% bifurcation with 'ode_HB2po'.
run_names.hopf_to_PO = 'run05_hopf_to_PO';

% Run continuation
hopf_to_PO;

%------------------------------------------%
%%     Compute Initial Periodic Orbit     %%
%------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.
run_names.initial_periodic_orbit = 'run06_initial_periodic_orbit';

% Run continuation
initial_periodic_orbit;

%-------------------------------------------------------------------------%
%%               Compute Floquet Bundle at Zero Phase Point              %%
%-------------------------------------------------------------------------%
% Here we compute the stable Floquet bundle of the periodic orbit, as well
% as the perpendicular vector, w.

% Add continuation scripts to path
addpath('./continuation_scripts/floquet_bundle/');

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Adjoint equations: Functions (for floquet_mu and floquet_wnorm)
% funcs.floquet = {@floquet_adjoint};
funcs.floquet = floquet_adjoint_symbolic();

% Boundary conditions: Floquet multipliers
% bcs_funcs.bcs_floquet = {@bcs_floquet};
bcs_funcs.bcs_floquet = bcs_floquet_symbolic();

%----------------------------------------------------------------%
%%     Compute Floquet Bundle at Zero Phase Point (with mu)     %%
%----------------------------------------------------------------%
% We now add the adjoint function and Floquet boundary conditions to
% compute the adjoint (left or right idk) eigenvectors and eigenvalues.
% This will give us the perpendicular vector to the tangent of the periodic
% orbit. However, this will only be for the eigenvector corresponding to
% the eigenvalue \mu = 1. Hence, here we continue in \mu (mu_f) until
% mu_f = 1.
run_names.compute_floquet_1 = 'run07_compute_floquet_bundle_1_mu';

% Run continuation script
floquet_mu;

%--------------------------------------------------------------------%
%%     Compute Floquet Bundle at Zero Phase Point (with w_norm)     %%
%--------------------------------------------------------------------%
% Having found the solution (branching point 'BP') corresponding to
% \mu = 1, we can continue in the norm of the vector w (w_norm), until the
% norm is equal to zero. Then we will have the correct perpendicular
% vector.
run_names.compute_floquet_2 = 'run08_compute_floquet_bundle_2_w';

% Run continuation script
floquet_wnorm;

%-------------------------------------------------------------------------%
%%                   Phase Response Curve Calculation                    %%
%-------------------------------------------------------------------------%
% We set up the phase resetting problem by creating four segments of the
% periodic orbit, with boundary conditions described in a paper somewhere.

% Add continuation scripts to path
addpath('./continuation_scripts/phase_reset/');
addpath('./plotting_scripts/phase_reset/');

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

% Boundary conditions: Segments 1 and 2
% bcs_funcs.bcs_seg1_seg2 = {@bcs_PR_seg1_seg2};
bcs_funcs.bcs_seg1_seg2 = bcs_PR_seg1_seg2_symbolic();

% Boundary conditions: Segment 3
% bcs_funcs.bcs_seg3 = {@bcs_PR_seg3};
bcs_funcs.bcs_seg3 = bcs_PR_seg3_symbolic();

% Boundary conditions: Segment 4
% bcs_funcs.bcs_seg4 = {@bcs_PR_seg4};
bcs_funcs.bcs_seg4 = bcs_PR_seg4_symbolic();

%------------------------------------------------------%
%%     First Continuation: Perturbation Amplitude     %%
%------------------------------------------------------%
% We compute the first phase resetting curve.
run_names.phase_reset_perturbation = 'run09_phase_reset_perturbation';

% Run continuation script
% phase_reset_1_perturbation

%--------------------------------------------------%
%%     First Continuation: Perturbation Angle     %%
%--------------------------------------------------%
% We compute the first phase resetting curve.
run_names.phase_reset_theta_perturb = 'run10_phase_reset_theta_perturb';

% Run continuation script
% phase_reset_2_theta_perturb
% phase_reset_2_theta_perturb_multi

%--------------------------------------------%
%%     Second Continuation: Compute PTCs     %
%--------------------------------------------%
% % Run name
% run_names.phase_transition_curve = 'run10_phase_reset_PTC_single';
% 
% % Run continuation script
% phase_reset_PTC;

% Run name
run_names.phase_transition_curve = 'run10_phase_reset_PTC_multi';

% Run continuation script
% phase_reset_PTC_scan;

%-------------------------------------------------------%
%%     First Continuation: theta_new and theta_old     %%
%-------------------------------------------------------%
% We compute the first phase resetting curve.
run_names.phase_reset_thetas = 'run11_phase_reset_thetas';

% Run continuation script
% follow_theta_new_theta_old;

%--------------------------------------------------------%
%%     Second Continuation: theta_new and theta_old     %%
%--------------------------------------------------------%
% We compute the first phase resetting curve.
run_names.isochron_test = 'run12_isochron_test';

% Run continuation script
% isochron_test;

%-------------------------------------------------------------------------%
%%                                PLOTS                                  %%
%-------------------------------------------------------------------------%
% Write data to txt files
% write_some_data_idk(run_names);
