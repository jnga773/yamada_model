%-------------------------------------------------------------------------%
%%                      Compute Singlular Amplitude                      %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.initial_periodic_orbit;
% Which run this continuation continues from
run_old = run_names.hopf_to_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');
label_old = label_old(1);

% label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
% label_old = max(label_old);

% Print to console
fprintf("~~~ Initial Periodic Orbit: Sixth Run (initial_periodic_orbit.m) ~~~ \n");
fprintf('Find new periodic orbit \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_soln = calculate_periodic_orbit(run_old, label_old);

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
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set initial guess to 'coll'
prob = ode_isol2coll(prob, 'initial_PO', funcs.field{:}, ...
                     data_soln.t, data_soln.x, pnames, data_soln.p);

% Add equilibrium points for non trivial steady states
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
prob = ode_ep2ep(prob, 'x0',   run_old, label_old);

% Glue parameters and apply boundary condition
prob = glue_parameters(prob, bcs_funcs.bcs_PO);

% Event for Nonlinear Photonics abstract
% prob = coco_add_event(prob, 'PO_PT', 'A', 7.3757);

% Event for A = 7.5
prob = coco_add_event(prob, 'PO_PT', 'A', data_soln.p(2));

% Run COCO
coco(prob, run_new, [], 1, {'A', 'gamma'}, A_range);

%-------------------------------------------------------------------------%
%%                             TESTING PLOTS                             %%
%-------------------------------------------------------------------------%
% Label for solution plot
label_plot = coco_bd_labs(coco_bd_read(run_new), 'PO_PT');
label_plot = label_plot(1);
% label_plot = 1;

% Calculate stable manifold of saddle point 'q' and save data to .mat in 
% ./data/ directory
save_initial_PO_data(run_new, label_plot);

% Plot solution
% plot_test_PO(run_new, label_plot);
plot_initial_periodic_orbit_COLL(save_figure);
plot_initial_periodic_orbit_2D_quasi(save_figure);
plot_initial_periodic_orbit_2D(save_figure);

%-------------------------------------------------------------------------%
%%                               FUNCTIONS                               %%
%-------------------------------------------------------------------------%
function data_out = calculate_periodic_orbit(run_in, label_in)
  % data_out = calculate_periodic_orbit(run_in, label_in)
  %
  % Calculates initial periodic solution using ode45

  %------------------%
  %     Read Data    %
  %------------------%
  % Read previous solution
  [sol, data] = po_read_solution('', run_in, label_in);

  % Initial period
  T_sol    = sol.T;
  % Parameters
  p_sol    = sol.p;
  % Read time data
  tbp      = sol.tbp;
  % Read state space solution
  xbp_read = sol.xbp;

  %----------------------------%
  %     Calculate Solution     %
  %----------------------------%
  % Need to find the point where e1 . F(x(t)) = 0, that is
  % the maximum of the first component.
  [~, max_idx] = max(xbp_read(:, 1));

  % Shift everything around
  if max_idx > 1
    x0 = [xbp_read(max_idx:end, :); xbp_read(2:max_idx, :)];
    t0 = [tbp(max_idx:end) - tbp(max_idx); tbp(2:max_idx) + (tbp(end) - tbp(max_idx))];
  end

  %----------------%
  %     Output     %
  %----------------%
  data_out.p      = p_sol;
  data_out.t      = t0;
  data_out.x      = x0;

end

function prob_out = glue_parameters(prob_in, bcs_PO_in)
  % prob_out = glue_parameter(prob_in)
  %
  % Glue the parameters of the EP segments and PO segment together 
  % (as they're all the same anyway)

  %---------------%
  %     Input     %
  %---------------%
  % Input continuation problem structure
  prob = prob_in;

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read index data for periodic orbit segment
  [data, uidx] = coco_get_func_data(prob, 'initial_PO.coll', 'data', 'uidx');

  % Read index data for equilibrium points
  [data1, uidx1] = coco_get_func_data(prob, 'xpos.ep', 'data', 'uidx');
  [data2, uidx2] = coco_get_func_data(prob, 'xneg.ep', 'data', 'uidx');
  [data3, uidx3] = coco_get_func_data(prob, 'x0.ep',   'data', 'uidx');

  % Index mapping
  maps     = data.coll_seg.maps;
  maps1    = data1.ep_eqn;
  maps2    = data2.ep_eqn;
  maps3    = data3.ep_eqn;

  %-----------------------------%
  %     Boundary Conditions     %
  %-----------------------------%
  % Apply periodic orbit boundary conditions and special phase condition
  prob = coco_add_func(prob, 'bcs_po', bcs_PO_in{:}, data, 'zero', 'uidx', ...
                       uidx([maps.x0_idx(1:data.xdim); ...
                             maps.x1_idx(1:data.xdim); ...
                             maps.p_idx(1:data.pdim)]));

  %-------------------------%
  %     Glue Parameters     %
  %-------------------------%
  prob = coco_add_glue(prob, 'glue_p1', uidx(maps.p_idx), uidx1(maps1.p_idx));
  prob = coco_add_glue(prob, 'glue_p2', uidx(maps.p_idx), uidx2(maps2.p_idx));
  prob = coco_add_glue(prob, 'glue_p3', uidx(maps.p_idx), uidx3(maps3.p_idx));

  %----------------%
  %     Output     %
  %----------------%
  prob_out = prob;

end

function save_initial_PO_data(run_in, label_in)
  % save_initial_PO_data(run_in, label_in)
  %
  % Reads periodic orbit solution data from COCO solution, calculates the
  % one-dimensional stable manifold of the "central" saddle point 'q', and
  % saves the data to './data/initial_PO.mat'.

  % Data matrix filename
  filename_out = './data/initial_PO.mat';

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Periodic orbit solution
  sol_PO = coll_read_solution('initial_PO', run_in, label_in);
  xbp_PO = sol_PO.xbp;
  p_PO   = sol_PO.p;

  % Equilibrium point solutions
  sol_pos = ep_read_solution('xpos', run_in, label_in);
  sol_neg = ep_read_solution('xneg', run_in, label_in);
  sol_0   = ep_read_solution('x0', run_in, label_in);

  x_pos = sol_pos.x;
  x_neg = sol_neg.x;
  x_0   = sol_0.x;

  %-----------------------------------%
  %     Calculate Stable Manifold     %
  %-----------------------------------%
  W_q_stable = calc_stable_manifold(run_in, label_in);

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_out, 'x_pos', 'x_neg', 'x_0', 'xbp_PO', 'p_PO', 'W_q_stable');

end

function x_out = calc_stable_manifold(run_in, label_in)
  % x_out = calc_stable_manifold(run_in, label_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read EP solution
  sol_EP   = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium point
  x_pos = sol_EP.x;
  % Parameters
  p     = sol_EP.p;
  
  %------------------------------%
  %     Calculate EigenStuff     %
  %------------------------------%
  % Jacobian
  J_stable = yamada_DFDX(x_pos, p);

  % Calculate eigenvalues and eigenvectors
  [eigvec, eigval] = eig(J_stable);

  % Indices for stable eigenvectors (eigval < 0)
  % stable_index = find(diag(eigval) < 0);
  stable_index = 3;

  % Stable eigenvector
  vec_s = eigvec(:, stable_index);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps1 = -0.01;
  % Time span
  % t_span1 = [0.0, -17.0];
  t_span1 = -16.5 : 0.01 : 0.0;
  t_span1 = flip(t_span1);

  % Initial vector
  x_init1 = x_pos + (eps1 * vec_s);

  % Integrate using ode45
  [~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span1, x_init1);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps2 = 0.01;
  % Time span
  % t_span2 = [0.0, -25.0];
  t_span2 = -25.0 : 0.01 : 0.0;
  t_span2 = flip(t_span2);

  % Initial vector
  x_init2 = x_pos + (eps2 * vec_s);

  % Integrate using ode45
  [~, W2] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span2, x_init2);

  %----------------%
  %     Output     %
  %----------------%
  x_out = [flip(W2); W1];

end
