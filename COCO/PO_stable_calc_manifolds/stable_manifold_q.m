%=========================================================================%
%                   YAMADA MODEL (Stable Manifold of q)                   %
%=========================================================================%
% We compute the stable manifold of the stable periodic orbit in Region 7,
% where there exists a stable/attracting periodic orbit. The Yamada mode
% equations are
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
addpath('./continuation_scripts/stable_manifold_q/');

% Add plotting scripts
addpath('./plotting_scripts/stable_manifold_q/');

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
bcs_funcs.bcs_initial = {@bcs_Wq_initial};

% Boundary conditions: Final condition
bcs_funcs.bcs_final = {@bcs_Wq_final};

%=========================================================================%
%                     CALCULATE STABLE MANIFOLD OF Q                      %
%=========================================================================%
% We compute the stable manifold of the saddle q in forward and backwards
% directions, based on original solution data from
% ./data_mat/initial_PO.mat'.
%-------------------------------------------------------------------------%
%%               Calculate Stable Manifold (q) Segments - 1              %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.stable_manifold1 = 'run06_stable_manifold_seg1';
run_new = run_names.stable_manifold1;

% Print to console
fprintf("~~~ Stable Manifold: Second Run ~~~ \n");
fprintf('Calculate one of the stable-manifold branches \n');
fprintf('Run name: %s \n', run_new);

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_isol = calc_initial_solution_Wsq('./data_mat/initial_PO.mat');

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

% Set step size
prob = coco_set(prob, 'cont', 'h_min', 5e-1, 'h', 1, 'h_max', 1e1);

% Set PtMX steps
PtMX = 50;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Continue periodic orbits
prob = ode_isol2coll(prob, 'initial_PO', funcs.field{:}, ...
                     data_isol.tbp_PO, data_isol.xbp_PO, data_isol.pnames, data_isol.p);

% Add collocation trajectory segment for stable manifold
prob = ode_isol2coll(prob, 'W1', funcs.field{:}, ...
                     data_isol.t0, data_isol.x_init_1, data_isol.p);
prob = ode_isol2coll(prob, 'W2', funcs.field{:}, ...
                     data_isol.t0, data_isol.x_init_2, data_isol.p);

% Continue equilibrium points for non trivial steady states
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   data_isol.xpos, data_isol.p);
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, ...
                   data_isol.xneg, data_isol.p);
prob = ode_isol2ep(prob, 'x0', funcs.field{:}, ...
                   data_isol.x0, data_isol.p);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, bcs_funcs, data_isol.eps);

%------------------------%
%     Add COCO Event     %
%------------------------%
prob = coco_add_event(prob, 'Del1', 'W_seg1', 0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[0.0, 10.0], [-1, 1e3], []};
coco(prob, run_new, [], 1, {'W_seg1', 'T1', 'W_seg2'}, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'Del1');

% Plot solution
plot_orbit_and_Wq_solution(run_new, label_plot);

%-------------------------------------------------------------------------%
%%               Calculate Stable Manifold (q) Segments - 2              %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.stable_manifold2 = 'run07_stable_manifold_seg2';
run_new = run_names.stable_manifold2;
run_old = run_names.stable_manifold1;

% Previous solution label
label_old = coco_bd_labs(coco_bd_read(run_old), 'Del1');

% Print to console
fprintf("~~~ Stable Manifold: Second Run ~~~ \n");
fprintf('Calculate one of the stable-manifold branches \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

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

% Set step size
prob = coco_set(prob, 'cont', 'h_min', 5e-1, 'h', 1, 'h_max', 1e1);

% Set PtMX steps
PtMX = 50;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Add collocation trajectory segment for stable manifold
prob = ode_coll2coll(prob, 'W1', run_old, label_old);
prob = ode_coll2coll(prob, 'W2', run_old, label_old);

% Continue periodic orbits
prob = ode_coll2coll(prob, 'initial_PO', run_old, label_old);

% Continue equilibrium points for non trivial steady states
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
prob = ode_ep2ep(prob, 'x0', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Grab the epsilong value
[data, chart] = coco_read_solution('bcs_final', run_old, label_old);
eps = chart.x(data.eps_idx);

% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, bcs_funcs, eps);

%------------------------%
%     Add COCO Event     %
%------------------------%
prob = coco_add_event(prob, 'Del2', 'W_seg2', 0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[-50, 0], [0.0, 1e3], []};
coco(prob, run_new, [], 1, {'W_seg2', 'T2', 'W_seg1'}, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'Del2');

% Plot solution
plot_orbit_and_Wq_solution(run_new, label_plot);

%-------------------------------------------------------------------------%
%%                Calculate Stable Manifold (q)- Close eps               %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.close_eps = 'run08_stable_manifold_close_eps';
run_new = run_names.close_eps;
run_old = run_names.stable_manifold2;

% Previous solution label
label_old = coco_bd_labs(coco_bd_read(run_old), 'Del2') ;

% Print to console
fprintf("~~~ Stable Manifold: Third Run ~~~ \n");
fprintf('Close the initial distance eps \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

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

% Set step size
prob = coco_set(prob, 'cont', 'h_min', 1e-2, 'h', 1e-1, 'h_max', 5e-1);

% Set PtMX steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Add collocation trajectory segment for stable manifold
prob = ode_coll2coll(prob, 'W1', run_old, label_old);
prob = ode_coll2coll(prob, 'W2', run_old, label_old);

% Continue periodic orbits
prob = ode_coll2coll(prob, 'initial_PO', run_old, label_old);

% Continue equilibrium points for non trivial steady states
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
prob = ode_ep2ep(prob, 'x0', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Grab the epsilong value
[data, chart] = coco_read_solution('bcs_final', run_old, label_old);
eps = chart.x(data.eps_idx);

% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_Wsq(prob, bcs_funcs, eps);

%------------------------%
%     Add COCO Event     %
%------------------------%
prob = coco_add_event(prob, 'EPS0', 'eps', [1e-8, 1e-9, 4.5e-10]);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[1e-8, eps], [], []};
coco(prob, run_new, [], 1, {'eps', 'T1', 'T2'}, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
% label_plot = coco_bd_labs(coco_bd_read(run_new), '');
% label_plot = max(label_plot) - 1;
label_plot = coco_bd_labs(coco_bd_read(run_new), 'EP');
label_plot = max(label_plot);

% Plot solution
plot_orbit_and_Wq_solution(run_new, label_plot);

%-------------------%
%     Save Data     %
%-------------------%
% Save solution
save_data_Wsq(run_new, label_plot, './data_mat/stable_manifold_q.mat');

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%