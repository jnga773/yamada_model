%=========================================================================%
%                YAMADA MODEL (STABLE AND UNSTABLE ORBITS)                %
%=========================================================================%
% We compute stable and unstable periodic orbits in region 8 of the
% two-parameter bifurcation diagram of the Yamada model:
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
% addpath('./functions/fields/hardcoded/');
addpath('./functions/fields/');
% Add boundary condition functions to path
% addpath('./functions/bcs/hardcoded/');
addpath('./functions/bcs/');
% Add SymCOCO files to path
addpath('./functions/symcoco/');

% Add continuations script functions to path
addpath('./continuation_scripts/initial_PO/');

% Add plotting scripts
addpath('./plotting_scripts/');

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

%-----------------------%
%     Problem Setup     %
%-----------------------%
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
% Yamada model equations
% funcs.field = {@yamada, @yamada_DFDX, @yamada_DFDP};
funcs.field = yamada_symbolic();

% Boundary conditions: Period
% bcs_funcs.bcs_T = {@bcs_T};
bcs_funcs.bcs_T = bcs_T_symbolic();

% Boundary conditions: Periodic orbit
% bcs_funcs.bcs_PO = {@bcs_PO};
bcs_funcs.bcs_PO = bcs_PO_symbolic();

% Boundary conditions: Eigenvalues and eigenvectors
bcs_funcs.bcs_eig = {@bcs_eig};

% Boundary conditions: Initial condition
bcs_funcs.bcs_initial = {@bcs_initial};

% Boundary conditions: Final condition
bcs_funcs.bcs_final = {@bcs_final};

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
%%                   Compute Initial Equilibrium Point                   %%
%-------------------------------------------------------------------------%
% We compute and continue the equilibrium point of the model using
% the 'EP' toolbox constructor 'ode_isol2ep'.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.initial_EP = 'run01_initial_EP';
run_new = run_names.initial_EP;

% Print to console
fprintf('~~~ Initial Periodic Orbit: First Run ~~~\n');
fprintf('Solve for initial solution of the equilibrium point\n')
fprintf('Run name: %s\n', run_new);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up COCO problem
prob = coco_prob();

% Set upper bound of continuation steps in each direction along solution
prob = coco_set(prob, 'cont', 'PtMX', 50);

% Set up isol2ep problem
prob = ode_isol2ep(prob, '', funcs.field{:}, x0, pnames, p0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, 'A', A_range);

%-------------------------------------------------------------------------%
%%                   Continuation from Branching Point                   %%
%-------------------------------------------------------------------------%
% We switch branches at the BP point to find the Hopf bifurcations.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.branching_point = 'run02_branching_point_continuation';
run_new = run_names.branching_point;
% Previous run name
run_old = run_names.initial_EP;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'BP');
label_old = label_old(1);

% Print to console
fprintf('~~~ Initial Periodic Orbit: Second run ~~~\n');
fprintf('Continue bifurcations from the branching point\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up COCO problem
prob = coco_prob();

% Set NAdapt to 1?
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set upper bound of continuation steps in each direction along solution
PtMX = 50;
prob = coco_set(prob, 'cont', 'PtMX', [PtMX, 0]);

% Continue from branching point
prob = ode_BP2ep(prob, '', run_old, label_old);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, 'A', A_range);

%-------------------------------------------------------------------------%
%%                           Move Hopf A Value                           %%
%-------------------------------------------------------------------------%
% Continuing from a Hopf bifurcation with 'ode_HB2HB', we vary
% the 'A' parameter to A = 7.3757

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.move_hopf = 'run03_move_hopf';
run_new = run_names.move_hopf;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'HB');
label_old = label_old(1);

% Print to console
fprintf("~~~ Initial Periodic Orbit: Third Run ~~~ \n");
fprintf('Move the gamma value \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 5e-2, 'h0', 5e-2, 'h_max', 5e-2);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set upper bound of continuation steps in each direction along solution
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Initial solution to periodic orbit (COLL Toolbox)
prob = ode_HB2HB(prob, '', run_old, label_old);

% Saved-point solution for A = 7.3757
prob = coco_add_event(prob, 'H_PT', 'A', 7.3757);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'}, A_range);

%-------------------------------------------------------------------------%
%%                        Hopf to Periodic Orbit                         %%
%-------------------------------------------------------------------------%
% Continue a family of periodic orbits emanating from the Hopf
% bifurcation with 'ode_HB2po'.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.hopf_to_PO = 'run05_hopf_to_PO';
run_new = run_names.hopf_to_PO;
% Which run this continuation continues from
run_old = run_names.move_hopf;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'H_PT');

% Print to console
fprintf("~~~ Initial Periodic Orbit: Fifth Run ~~~ \n");
fprintf('Periodic orbits from Hopf bifurcation \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------%
%     Calculate Things     %
%--------------------------%
% Read previous solution
sol = ep_read_solution('', run_old, label_old);

% Calculate non-trivial solutions
[xpos, xneg] = non_trivial_ss(sol.p);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 25);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% % Set step sizes
% h_size = 1e0;
% prob = coco_set(prob, 'cont', 'h_min', h_size, 'h0', h_size, 'h_max', h_size);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Continue from Hopf bifurcation
prob = ode_HB2po(prob, '', run_old, label_old);

% Follow non trivial solutions
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   xpos, sol.p);
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, ...
                   xneg, sol.p);
prob = ode_isol2ep(prob, 'x0',   funcs.field{:}, ...
                   x0,   sol.p);

% Glue parameters (defined in './continuation_scripts/glue_parameters.m')
prob = glue_parameters_PO(prob);

% Saved point for solution for gamma = 3.54e-2
prob = coco_add_event(prob, 'PO_PT', 'gamma', 3.54e-2);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
bd = coco(prob, run_new, [], 1, {'gamma', 'A'}, gamma_range);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution to plot
label_plot = sort(coco_bd_labs(coco_bd_read(run_new), 'PO_PT'));

% Create plots
plot_hopf_to_PO_solution(run_new, label_plot);

%-------------------------------------------------------------------------%
%%                   Re-Solve for Rotated Perioid Orbit                  %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.initial_PO = 'run06_initial_periodic_orbit';
run_new = run_names.initial_PO;
% Which run this continuation continues from
run_old = run_names.hopf_to_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');

% Print to console
fprintf("~~~ Initial Periodic Orbit: Sixth Run ~~~ \n");
fprintf('Find new periodic orbit \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_stable = calc_PO_initial_solution(run_old, label_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdpat
% prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Stable periodic orbit
prob = ode_isol2coll(prob, 'PO_stable', funcs.field{:}, ...
                     data_stable.t, data_stable.x, data_stable.pnames, data_stable.p, ...
                     '-var', eye(3));

% Add equilibrium points for non trivial steady states
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   data_stable.xpos, data_stable.p);
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, ...
                   data_stable.xneg, data_stable.p);
prob = ode_isol2ep(prob, 'x0', funcs.field{:}, ...
                   data_stable.x0, data_stable.p);

% Glue parameters and apply boundary condition
prob = apply_PO_boundary_conditions(prob, bcs_funcs.bcs_PO);

% Event for A = 7.5
prob = coco_add_event(prob, 'PO_PT', 'A', data_stable.p(2));

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'}, A_range);

%----------------------%
%    Testing Plots     %
%----------------------%
% Label for solution plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'PO_PT');
label_plot = label_plot(1);

% Plot solution
plot_initial_periodic_orbit_COLL(run_new, label_plot);

% Save solution
save_initial_PO_data(run_new, label_plot);

%--------------------------%
%     Monodromy Matrix     %
%--------------------------%
label_plot = coco_bd_labs(coco_bd_read(run_new), 'PO_PT');

[sol, ~] = coll_read_solution('PO_stable', run_new, label_plot);
xbp_PO = sol.xbp;
fprintf('Head-point = (%.3f, %.3f, %.3f)\n', xbp_PO(1, :));

% Monodromy matrix
chart = coco_read_solution('PO_stable.coll.var', run_new, label_plot, 'chart');
data  = coco_read_solution('PO_stable.coll', run_new, label_plot, 'data');

% Create monodrony matrix
M1 = chart.x(data.coll_var.v1_idx);

fprintf('~~~ Monodromy Matrix ~~~\n');
fprintf('(%.7f, %.7f, %.7f)\n', M1(1, :));
fprintf('(%.7f, %.7f, %.7f)\n', M1(2, :));
fprintf('(%.7f, %.7f, %.7f)\n\n', M1(3, :));

% Get eigenvalues and eigenvectors of the Monodromy matrix
[vec_floquet, eig_floquet] = eig(M1);

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%