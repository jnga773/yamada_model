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
addpath('./continuation_scripts/initial_PO/');
addpath('./continuation_scripts/phase_reset/');
addpath('./continuation_scripts/stable_manifold/');

% Add plotting scripts
addpath('./plotting_scripts/');

%--------------------%
%     Parameters     %
%--------------------%
% Because we will only be looking at the (A, \gamma) plane, we will be
% setting values for a and B.
B = 5.8;
a = 1.8;

% Parameters for the periodic orbit
gamma_PO = 3.5e-2;
A_PO     = 7.4;

%-----------------------%
%     Problem Setup     %
%-----------------------%
% Parameter names
pnames = {'gamma', 'A', 'B', 'a'};

% Initial parameter values
p0 = [gamma_PO; A_PO; B; a];

% Initial state values
x0 = [10; 10; 10];

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

% Phase Reset Segment 1: Functions
% func.seg1 = {@func_seg1};
funcs.seg1 = func_seg1_symbolic();

% Phase Reset: Segment 2
% funcs.seg2 = {@func_seg2};
funcs.seg2 = func_seg2_symbolic();

% Boundary conditions: Phase-resetting segments
% bcs_funcs.bcs_PR = {@bcs_PR};
bcs_funcs.bcs_PR = bcs_PR_symbolic();

% Boundary conditions: W^{s}(q) initial condition
bcs_funcs.bcs_Wsq_initial = {@bcs_Wsq_initial};

% Boundary conditions: W^{s}(q) final condition
bcs_funcs.bcs_Wsq_final = {@bcs_Wsq_final};

%=========================================================================%
%%                   CALCULATE INITIAL PERIODIC ORBIT                    %%
%=========================================================================%
% We compute a family of periodic orbits, emanating from a Hopf
% bifurcation. We first compute the Hopf bifurcation line using the 'EP'
% toolbox, and then a family of periodic orbits with the 'PO' toolbox.
% Finally, we shift the state-space solution data such that, at t=0,
%                       x1(t=0) = max(x1) .
% We then verify this solution using the 'COLL' toolbox.

%-------------------------------------------------------------------------%
%%                 Confirm ODE45 Periodic Orbit Solution                 %%
%-------------------------------------------------------------------------%
% Calculate the periodic orbit using MATLAB's ode45 function.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.ode45_PO = 'run01_initial_PO_ode45';
run_new = run_names.ode45_PO;

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Initial Periodic Orbit: First Run\n');
fprintf(' Find new periodic orbit\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Continuation parameters : %s\n', 'A, gamma');
fprintf(' =====================================================================\n');

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_ode45 = calc_initial_solution_ODE45(x0, p0, funcs.field);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set initial guess to 'coll'
prob = ode_isol2po(prob, '', funcs.field{:}, ...
                   data_ode45.t, data_ode45.x, pnames, p0);

% Add equilibrium points for non trivial steady states
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   data_ode45.xpos, p0);
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, ...
                   data_ode45.xneg, p0);
prob = ode_isol2ep(prob, 'x0', funcs.field{:}, ...
                   data_ode45.x0, p0);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Glue parameters and apply boundary condition
prob = glue_parameters_PO(prob);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
prob = coco_add_event(prob, 'PO_PT', 'A', A_PO);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
bdtest = coco(prob, run_new, [], 1, {'A', 'gamma'});

%-------------------------------------------------------------------------%
%%                   Re-Solve for Rotated Perioid Orbit                  %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.initial_PO = 'run02_initial_periodic_orbit';
run_new = run_names.initial_PO;
% Which run this continuation continues from
run_old = run_names.ode45_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Initial Periodic Orbit: Second Run\n');
fprintf(' Rotate periodic orbit\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, gamma');
fprintf(' =====================================================================\n');

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_PO = calc_initial_solution_PO(run_old, label_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set initial guess to 'coll'
prob = ode_isol2coll(prob, 'initial_PO', funcs.field{:}, ...
                     data_PO.t, data_PO.x, pnames, data_PO.p);

% Add equilibrium points for non trivial steady states
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
prob = ode_ep2ep(prob, 'x0',   run_old, label_old);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_PO(prob, bcs_funcs.bcs_PO);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Event for A = 7.5
prob = coco_add_event(prob, 'PO_PT', 'A', data_PO.p(2));

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'});

%-------------------%
%     Save Data     %
%-------------------%
% Label of solution to save
label_plot = coco_bd_labs(coco_bd_read(run_new), 'PO_PT');

% Plot solution
plot_initial_PO(run_new, label_plot);

%=========================================================================%
%%                    CALCULATE PHASE RESET SOLUTIONS                    %%
%=========================================================================%
% We set up the phase resetting segments, but for a zero-amplitude
% perturbation. We then continue along the two phase to find points along
% the periodic orbit which pass the equilibrium point q in the X-Z plane.

%-------------------------------------------------------------------------%
%%                        Move Along Period Orbit                        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.PR_move_theta = 'run03_PR_move_theta';
run_new = run_names.PR_move_theta;
% Which run this continuation continues from
run_old = run_names.initial_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Phase Transition Curve: First Run\n');
fprintf(' Move theta along Gamma\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'theta_gamma');
fprintf(' =====================================================================\n');

%-------------------%
%     Read Data     %
%-------------------%
% Set initial conditions from previous solutions
data_PR = calc_initial_solution_PR(run_old, label_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set number of steps
PtMX = 1000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

% Set norm to int
prob = coco_set(prob, 'cont', 'norm', inf);

%------------------%
%     Set NTST     %
%------------------%
% Set NTST for the orbit segments
NTST = 50;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST);
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST);

%------------------------------------%
%     Create Trajectory Segments     %
%------------------------------------%
% Segment 1
prob = ode_isol2coll(prob, 'seg1', funcs.seg1{:}, ...
                     data_PR.t_seg1, data_PR.x_seg1, data_PR.pnames, data_PR.p0);
% Segment 2
prob = ode_isol2coll(prob, 'seg2', funcs.seg2{:}, ...
                     data_PR.t_seg2, data_PR.x_seg2, data_PR.p0);

% Add equilibrium point
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   data_PR.xpos, data_PR.p0(1:data_PR.pdim));

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply all boundary conditions, glue parameters together, and
% all that other good COCO stuff. Looking the function file
% if you need to know more ;)
prob = apply_boundary_conditions_PR(prob, data_PR, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% We monitor the distance of the second component of \gamma_{\theta_{old}}
% from the equilibrium point q, which is the point we want to pass.
prob = coco_add_event(prob, 'q_PT', 'q_dist', 0.0);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and parameter range
pcont = {'theta_gamma', 'q_dist'};
prange = {[0.0, 1.0], []};

% Run COCO
coco(prob, run_new, [], 1, pcont, prange);

%=========================================================================%
%%                    CALCULATE STABLE MANIFOLD OF Q                     %%
%=========================================================================%
% We compute the stable manifold of the saddle q in forward and backwards
% directions.

%-------------------------------------------------------------------------%
%%               Calculate Stable Manifold (Part 1) - Run 1              %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.grow_Wsq_part1 = {'run04_grow_Wsq', 'part1'};
run_new = run_names.grow_Wsq_part1;
% Which run this continuation continues from
run_old = run_names.PR_move_theta;

% Continuation point
label_old = 1;

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Stable Manifold of q: First Run\n');
fprintf(' Calculate one of the stable-manifold branches \n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : {%s, %s}\n', run_new{1}, run_new{2});
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'Wsq_dist, T_Wsq, theta_gamma');
fprintf(' =====================================================================\n');

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Initial spacing
epsilon = -1e-3;

% Calculate dem tings
data_Wsq = calc_initial_solution_Wsq(run_old, label_old, epsilon, funcs);
data_Wsq.p_maps = data_PR.p_maps;

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2);
prob = coco_set(prob, 'cont', 'h0', 1e-1);
prob = coco_set(prob, 'cont', 'h_max', 1e1);

% Set number of steps
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

%------------------%
%     Set NTST     %
%------------------%
% Set NTST for the orbit segments
NTST = 50;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST);
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST);
prob = coco_set(prob, 'Wsq.coll', 'NTST', NTST);

%--------------------------------------%
%     Continue Trajectory Segments     %
%--------------------------------------%
% Segment 1
prob = ode_coll2coll(prob, 'seg1', run_old, label_old);
% Segment 2
prob = ode_coll2coll(prob, 'seg2', run_old, label_old);

% Add equilibrium point
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-------------------------------------------%
%     Setup Manifold Trajectory Segment     %
%-------------------------------------------%
% Set initial guess for the manifold segment
prob = ode_isol2coll(prob, 'Wsq', funcs.field{:}, ...
                     data_Wsq.t0, data_Wsq.x0, data_Wsq.p0(1:data_Wsq.pdim));

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, data_Wsq, bcs_funcs, data_Wsq.epsilon);

%------------------------%
%     Add CoCo Event     %
%------------------------%
prob = coco_add_event(prob, 'W_HIT', 'Wsq_dist', 0);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and parameter range
pcont  = {'Wsq_dist', 'T_Wsq', 'theta_gamma'};
prange = {[-100, 0], [0, 1e6], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'W_HIT');

% Plot solution
plot_orbit_and_Wq_solution(run_new, label_plot);

%-------------------------------------------------------------------------%
%%               Calculate Stable Manifold (Part 1) - Run 2              %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.follow_theta_Wsq_part1 = {'run05_follow_theta_Wsq', 'part1'};
run_new = run_names.follow_theta_Wsq_part1;
% Which run this continuation continues from
run_old = run_names.grow_Wsq_part1;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'W_HIT');

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Stable Manifold of q: Third Run\n');
fprintf(' Move theta_old along Gamma\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : {%s, %s}\n', run_new{1}, run_new{2});
fprintf(' Previous run name       : {%s, %s}\n', run_old{1}, run_old{2});
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'theta_gamma, T_Wsq, A_perturb, theta_perturb');
fprintf(' =====================================================================\n');

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-4);
prob = coco_set(prob, 'cont', 'h0', 1e-3);
prob = coco_set(prob, 'cont', 'h_max', 1e0);

% Set number of steps
PtMX = 2000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

%------------------%
%     Set NTST     %
%------------------%
% Set NTST for the orbit segments
NTST = 50;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST);
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST);
prob = coco_set(prob, 'Wsq.coll', 'NTST', NTST);

%--------------------------------------%
%     Continue Trajectory Segments     %
%--------------------------------------%
% Segment 1
prob = ode_coll2coll(prob, 'seg1', run_old, label_old);
% Segment 2
prob = ode_coll2coll(prob, 'seg2', run_old, label_old);

% Add equilibrium point
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-------------------------------------------%
%     Setup Manifold Trajectory Segment     %
%-------------------------------------------%
% Set initial guess for the manifold segment
prob = coco_set(prob, 'Wsq.coll', 'NTST', 100);
prob = ode_coll2coll(prob, 'Wsq', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Grab the epsilong value
[data_bcs, chart_bcs] = coco_read_solution('bcs_Wsq_final', run_old, label_old);
epsilon_read = chart_bcs.x(data_bcs.eps_idx);

% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, data_Wsq, bcs_funcs, epsilon_read);

%------------------------%
%     Add CoCo Event     %
%------------------------%
% Add saved point for when theta_old crosses over q
prob = coco_add_event(prob, 'q_PT', 'q_dist', 0.0);
% Add saved points for theta_perturb = 0, 45, and 90 degrees
SP_points = [0.0, 0.125, 0.25, 1.0, 1.125, 1.25];
prob = coco_add_event(prob, 'SP', 'theta_perturb', SP_points);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and parameter range
pcont  = {'theta_gamma', 'T_Wsq', 'A_perturb', 'theta_perturb'};
prange = {[0.0, 2.0], [0, 1e6], [], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% % Solution label to plot
% label_plot = 4;
% 
% % Plot solution
% plot_orbit_and_Wq_solution(run_new, label_plot);

% Plot A_perturb
% plot_run_perturbation_two_axis(run_new);
plot_bif_run_A_perturb(run_new);
plot_bif_run_theta_old(run_new);

%-------------------------------------------------------------------------%
%%               Calculate Stable Manifold (Part 2) - Run 1              %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.grow_Wsq_part2 = {'run04_grow_Wsq', 'part2'};
run_new = run_names.grow_Wsq_part2;
% Which run this continuation continues from
run_old = run_names.PR_move_theta;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'q_PT');
label_old = min(label_old) + 1;

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Stable Manifold of q: First Run\n');
fprintf(' Calculate one of the stable-manifold branches \n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : {%s, %s}\n', run_new{1}, run_new{2});
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'Wsq_dist, T_Wsq, theta_gamma');
fprintf(' =====================================================================\n');

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Initial spacing
epsilon = 1e-3;

% Calculate dem tings
data_Wsq = calc_initial_solution_Wsq(run_old, label_old, epsilon, funcs);
data_Wsq.p_maps = data_PR.p_maps;

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2);
prob = coco_set(prob, 'cont', 'h0', 1e-1);
prob = coco_set(prob, 'cont', 'h_max', 1e1);

% Set number of steps
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

%------------------%
%     Set NTST     %
%------------------%
% Set NTST for the orbit segments
NTST = 50;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST);
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST);
prob = coco_set(prob, 'Wsq.coll', 'NTST', NTST);

%--------------------------------------%
%     Continue Trajectory Segments     %
%--------------------------------------%
% Segment 1
prob = ode_coll2coll(prob, 'seg1', run_old, label_old);
% Segment 2
prob = ode_coll2coll(prob, 'seg2', run_old, label_old);

% Add equilibrium point
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-------------------------------------------%
%     Setup Manifold Trajectory Segment     %
%-------------------------------------------%
% Set initial guess for the manifold segment
prob = ode_isol2coll(prob, 'Wsq', funcs.field{:}, ...
                     data_Wsq.t0, data_Wsq.x0, data_Wsq.p0(1:data_Wsq.pdim));

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, data_Wsq, bcs_funcs, data_Wsq.epsilon);

%------------------------%
%     Add CoCo Event     %
%------------------------%
prob = coco_add_event(prob, 'W_HIT', 'Wsq_dist', 0);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and parameter range
pcont  = {'Wsq_dist', 'T_Wsq', 'theta_gamma'};
prange = {[0, 100], [0, 1e6], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'W_HIT');

% Plot solution
plot_orbit_and_Wq_solution(run_new, label_plot);


%-------------------------------------------------------------------------%
%%               Calculate Stable Manifold (Part 2) - Run 2              %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.follow_theta_Wsq_part2 = {'run05_follow_theta_Wsq', 'part2'};
run_new = run_names.follow_theta_Wsq_part2;
% Which run this continuation continues from
run_old = run_names.grow_Wsq_part2;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'W_HIT');

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Stable Manifold of q: Third Run\n');
fprintf(' Move theta_old along Gamma\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : {%s, %s}\n', run_new{1}, run_new{2});
fprintf(' Previous run name       : {%s, %s}\n', run_old{1}, run_old{2});
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'theta_gamma, T_Wsq, A_perturb, theta_perturb');
fprintf(' =====================================================================\n');

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-4);
prob = coco_set(prob, 'cont', 'h0', 1e-3);
prob = coco_set(prob, 'cont', 'h_max', 1e0);

% Set number of steps
PtMX = 2000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

%------------------%
%     Set NTST     %
%------------------%
% Set NTST for the orbit segments
NTST = 50;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST);
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST);
prob = coco_set(prob, 'Wsq.coll', 'NTST', NTST);

%--------------------------------------%
%     Continue Trajectory Segments     %
%--------------------------------------%
% Segment 1
prob = ode_coll2coll(prob, 'seg1', run_old, label_old);
% Segment 2
prob = ode_coll2coll(prob, 'seg2', run_old, label_old);

% Add equilibrium point
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-------------------------------------------%
%     Setup Manifold Trajectory Segment     %
%-------------------------------------------%
% Set initial guess for the manifold segment
prob = coco_set(prob, 'Wsq.coll', 'NTST', 100);
prob = ode_coll2coll(prob, 'Wsq', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Grab the epsilong value
[data_bcs, chart_bcs] = coco_read_solution('bcs_Wsq_final', run_old, label_old);
epsilon_read = chart_bcs.x(data_bcs.eps_idx);

% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, data_Wsq, bcs_funcs, epsilon_read);

%------------------------%
%     Add CoCo Event     %
%------------------------%
% Add saved point for when theta_old crosses over q
prob = coco_add_event(prob, 'q_PT', 'q_dist', 0.0);
% Add saved points for theta_perturb = 0, 45, and 90 degrees
SP_points = [0.0, 0.125, 0.25, 1.0, 1.125, 1.25];
prob = coco_add_event(prob, 'SP', 'theta_perturb', SP_points);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and parameter range
pcont  = {'theta_gamma', 'T_Wsq','A_perturb', 'theta_perturb'};
prange = {[0.0, 2.0], [0, 1e6], [], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% % Solution label to plot
% label_plot = 4;
% 
% % Plot solution
% plot_orbit_and_Wq_solution(run_new, label_plot);

% Plot A_perturb
% plot_run_perturbation_two_axis(run_new);
plot_bif_run_A_perturb(run_new);
plot_bif_run_theta_old(run_new);

%=========================================================================%
%%                          COMBINE DATA TO PLOT                         %%
%=========================================================================%
% Get run names
data_run1 = run_names.follow_theta_Wsq_part1;
data_run2 = run_names.follow_theta_Wsq_part2;
data_runq = run_names.PR_move_theta;

% Plot
plot_two_axis_figure(run_names.follow_theta_Wsq_part1, ...
                     run_names.follow_theta_Wsq_part2, ...
                     run_names.PR_move_theta);

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%
