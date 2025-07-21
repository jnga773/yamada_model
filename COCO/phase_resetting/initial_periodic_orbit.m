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
addpath('./continuation_scripts/initial_periodic_orbit/');

% Add plotting scripts
addpath('./plotting_scripts/initial_periodic_orbit');

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

% Adjoint equations: Functions (for floquet_mu and floquet_wnorm)
% funcs.VAR = {@VAR_symbolic};
funcs.VAR = VAR_symbolic();

% Boundary conditions: Floquet multipliers
% bcs_funcs.bcs_VAR = {@bcs_VAR};
bcs_funcs.bcs_VAR = bcs_VAR_symbolic();

%=========================================================================%
%                    CALCULATE INITIAL PERIODIC ORBIT                     %
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
run_names.initial_PO_ode45 = 'run01_initial_PO_ode45';
run_new = run_names.initial_PO_ode45;

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
run_names.initial_PO = 'run02_initial_PO';
run_new = run_names.initial_PO;
% Which run this continuation continues from
run_old = run_names.initial_PO_ode45;

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

% Calculate stable manifold of saddle point 'q' and save data to .mat in 
% ./data_mat/ directory
save_data_PO(run_new, label_plot, './data_mat/initial_PO.mat');

% Plot solution
plot_initial_PO(run_new, label_plot);

%=========================================================================%
%%            Compute Floquet Bundle at Zero Phase Point (mu)            %%
%=========================================================================%
% We now add the adjoint function and Floquet boundary conditions to
% compute the adjoint (left or right idk) eigenvectors and eigenvalues.
% This will give us the perpendicular vector to the tangent of the periodic
% orbit. However, this will only be for the eigenvector corresponding to
% the eigenvalue \mu = 1. Hence, here we continue in \mu (mu_s) until
% mu_s = 1.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.VAR_mu = 'run03_VAR_mu';
run_new = run_names.VAR_mu;
% Which run this continuation continues from
run_old = run_names.initial_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Floquet Bundle: First Run\n');
fprintf(' Calculate stable Floquet bundle eigenvalue\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'mu_s, w_norm');
fprintf(' =====================================================================\n');

%--------------------------%
%     Calculate Things     %
%--------------------------%
data_adjoint = calc_initial_solution_VAR(run_old, label_old);

%------------------------------------%
%     Setup Floquet Continuation     %
%------------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2, 'h0', 1e-2, 'h_max', 1e-2);

% Set PtMX
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set NTST
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdapt
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

% Add segment as initial solution
prob = ode_isol2coll(prob, 'adjoint', funcs.VAR{:}, ...
                     data_adjoint.t0, data_adjoint.x0, ...
                     data_adjoint.pnames, data_adjoint.p0);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply boundary conditions
prob = apply_boundary_conditions_VAR(prob, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Add event
prob = coco_add_event(prob, 'mu=1', 'mu_s', 1.0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'mu_s', 'w_norm'} , [0.0, 1.1]);

%-------------------------------------------------------------------------%
%%          Compute Floquet Bundle at Zero Phase Point (w_norm)          %%
%-------------------------------------------------------------------------%
% Having found the solution (branching point 'BP') corresponding to
% \mu = 1, we can continue in the norm of the vector w (w_norm), until the
% norm is equal to zero. Then we will have the correct perpendicular
% vector.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.VAR_wnorm = 'run04_VAR_wnorm';
run_new = run_names.VAR_wnorm;
% Which run this continuation continues from
run_old = run_names.VAR_mu;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'BP');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Floquet Bundle: Second Run\n');
fprintf(' Grow norm of stable Floquet bundle vector\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'mu_s, w_norm');
fprintf(' =====================================================================\n');

%------------------------------------%
%     Setup Floquet Continuation     %
%------------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set number of PtMX steps
PtMX = 2000;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set number of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Continue coll from previous branching point
prob = coco_set(prob, 'cont', 'branch', 'switch');
prob = ode_coll2coll(prob, 'adjoint', run_old, label_old);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply boundary conditions
prob = apply_boundary_conditions_VAR(prob, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Add event when w_norm = 1
prob = coco_add_event(prob, 'NORM1', 'w_norm', 1.0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'mu_s', 'w_norm'}, {[], [-1e-4, 1.1]});

%-------------------%
%     Save Data     %
%-------------------%
label_plot = coco_bd_labs(coco_bd_read(run_new), 'NORM1');

% Save solution to .mat to be read in 'yamada_PTC.m'
save_data_VAR(run_new, label_plot, './data_mat/solution_VAR.mat');

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%
