%=========================================================================%
%                YAMADA MODEL (Stable Manifold of Gamma_u)                %
%=========================================================================%
% We compute the stable manifold of the attracting periodic orbit in Region
% 7, where there exists a stable/attracting periodic orbit. The Yamada
% model equations are
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
addpath('./continuation_scripts/stable_manifold_PO/');

% Add plotting scripts
addpath('./plotting_scripts/stable_manifold_PO/');

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Yamada model equations
% funcs.field = {@yamada, @yamada_DFDX, @yamada_DFDP};
funcs.field   = yamada_symbolic();

% Boundary conditions: Periodic orbit
% bcs_funcs.bcs_PO = {@bcs_PO};
bcs_funcs.bcs_PO = bcs_PO_symbolic();

% Boundary conditions: Eigenvalues and eigenvectors
bcs_funcs.bcs_eig = {@bcs_eig};

% Boundary conditions: Initial condition
bcs_funcs.bcs_initial = {@bcs_W_PO_initial};

% Boundary conditions: Final condition
bcs_funcs.bcs_final = {@bcs_W_PO_final};

%=========================================================================%
%            CALCULATE STABLE MANIFOLD OF SADDLE PERIODIC ORBIT           %
%=========================================================================%
% We compute the stable manifold of the saddle q in forward and backwards
% directions, based on original solution data from
% ./data_mat/initial_PO.mat'.

%-------------------------------------------------------------------------%
%%                 Calculate Stable Manifold Segments - 1                %%
%-------------------------------------------------------------------------%
% Using previous parameters and sasved solutions, we compute the stable and
% saddle periodic orbits, the three stationary points, and two separate
% segments for the stable manifold of the periodic orbit.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.stable_manifold1 = 'run01_Ws_PO';
run_new = run_names.stable_manifold1;

% Print to console
fprintf("~~~ Stable Manifold (PO): First Run ~~~ \n");
fprintf('Calculate one of the stable-manifold branches \n');
fprintf('Run name: %s \n', run_new);

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_isol = calc_initial_solution_WsPO('./data_mat/initial_PO.mat');

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
% prob = coco_set(prob, 'cont', 'h_min', 1, 'h', 5, 'h_max', 1e1);

% Set PtMX steps
PtMX = 400;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Continue periodic orbits
prob = ode_isol2coll(prob, 'initial_PO', funcs.field{:}, ...
                     data_isol.tbp_PO, data_isol.xbp_PO, data_isol.pnames, data_isol.p, ...
                     '-var', eye(3));

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
prob = apply_boundary_conditions_WsPO(prob, bcs_funcs, data_isol.eps, data_isol.vec_s, data_isol.lam_s);

%------------------------%
%     Add CoCo Event     %
%------------------------%
prob = coco_add_event(prob, 'Del1', 'W_seg1', 0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[0, 20], [0.0, 1e4], []};
coco(prob, run_new, [], 1, {'W_seg1', 'T1', 'W_seg2'}, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'Del1');

% Plot solution
plot_orbit_and_W_PO_solution(run_new, label_plot);

%-------------------------------------------------------------------------%
%%                 Calculate Stable Manifold Segments - 2                %%
%-------------------------------------------------------------------------%
% Continuing from the previous solution where the W1 manifold segment was
% computed up to some boundary plane, we compute the W2 manifold segment.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.stable_manifold2 = 'run02_stable_manifold_seg2';
run_new = run_names.stable_manifold2;
run_old = run_names.stable_manifold1;

% Previous solution label
label_old = coco_bd_labs(coco_bd_read(run_old), 'Del1');

% Print to console
fprintf("~~~ Stable Manifold (q): Second Run ~~~ \n");
fprintf('Calculate one of the stable-manifold branches \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% % Set step size
% prob = coco_set(prob, 'cont', 'h_min', 5e-1, 'h', 1e0, 'h_max', 1e1);

% Set PtMX steps
PtMX = 300;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 25);

% Continue periodic orbits
prob = ode_coll2coll(prob, 'initial_PO', run_old, label_old, ...
                     '-var', eye(3));

% Add collocation trajectory segment for stable manifold
prob = ode_coll2coll(prob, 'W1', run_old, label_old);
prob = ode_coll2coll(prob, 'W2', run_old, label_old);

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

% Grab the Floquet eigenvector and eigenvalue
[data, chart] = coco_read_solution('bcs_eig', run_old, label_old);
vec_s = chart.x(data.vec_floquet_idx);
lam_s = chart.x(data.lam_floquet_idx);

% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_WsPO(prob, bcs_funcs, eps, vec_s, lam_s);

%------------------------%
%     Add CoCo Event     %
%------------------------%
prob = coco_add_event(prob, 'Del2', 'W_seg2', 0.0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[-100, 0], [0.0, 1e3], []};
coco(prob, run_new, [], 1, {'W_seg2', 'T2', 'T1'}, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'Del2');

% Plot solution
plot_orbit_and_W_PO_solution(run_new, label_plot);

%-------------------------------------------------------------------------%
%%                     Sweep Out From Periodic Orbit                     %%
%-------------------------------------------------------------------------%
% We now sweep through values of epsilon (eps), moving closer to and
% further away form the periodic orbit. 

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.stable_manifold3 = 'run03_stable_manifold_whanau';
run_new = run_names.stable_manifold3;
run_old = run_names.stable_manifold2;

% Previous solution label
label_old = coco_bd_labs(coco_bd_read(run_old), 'Del2');

% Print to console
fprintf("~~~ Stable Manifold (q): Third Run ~~~ \n");
fprintf('Calculate one of the stable-manifold branches \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 70);

% Set NAdpat
% prob = coco_set(prob, 'cont', 'NAdapt', 20);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set step size
% prob = coco_set(prob, 'cont', 'h_min', 1e-2, 'h', 5e-2, 'h_max', 1e-1);

% Set PtMX steps
PtMX = 2000;
prob = coco_set(prob, 'cont', 'PtMX', [PtMX, 0]);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Set norm
% prob = coco_set(prob, 'cont', 'norm', inf);

% Continue periodic orbits
prob = ode_coll2coll(prob, 'initial_PO', run_old, label_old, ...
                     '-var', eye(3));

% Add collocation trajectory segment for stable manifold
prob = ode_coll2coll(prob, 'W1', run_old, label_old);
prob = ode_coll2coll(prob, 'W2', run_old, label_old);

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

% Grab the Floquet eigenvector and eigenvalue
[data, chart] = coco_read_solution('bcs_eig', run_old, label_old);
vec_s = chart.x(data.vec_floquet_idx);
lam_s = chart.x(data.lam_floquet_idx);

% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_WsPO(prob, bcs_funcs, eps, vec_s, lam_s);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
prange = {[1e-10, eps], [], []};
coco(prob, run_new, [], 1, {'eps', 'T1', 'T2'}, prange);

%----------------------%
%    Testing Plots     %
%----------------------%
% Solution label to plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'EP');
label_plot = max(label_plot);
% label_plot = 12;

% Plot singlesolution
plot_orbit_and_W_PO_solution(run_new, label_plot);
% Plot solution scan
plot_solutions_scan(run_new)

%-------------------%
%     Save Data     %
%-------------------%
% Save solutions
save_data_WsPO(run_new, './data_mat/stable_manifold_PO.mat');

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%
