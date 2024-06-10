%-------------------------------------------------------------------------%
%%                     Parametrise the Heteroclinic                      %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.lins_method.continue_homoclinics;
% Which run this continuation continues from
run_old = run_names.lins_method.close_lingap;
% run_old = run_names.lins_method.close_eps2;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'Lin0');
% label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = max(label_old);

% Print to console
fprintf("~~~ Lin's Method: Sixth Run (ode_coll2coll) ~~~ \n");
fprintf('Continue constrained segments to find parametrisation of homoclinic \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------------------------%
%     Read x-Solutions from Previous Run     %
%--------------------------------------------%
single_segment = merge_segments(run_old, label_old);

% Extract stored deviations of heteroclinic trajectory end points from corresponding equilibria
[data, chart] = coco_read_solution('bcs_initial', run_old, label_old);
epsilon = chart.x(data.epsilon_idx);

% Eigenvalues and eigenvectors
[eigvecs, eigvals] = read_eigen_data(run_old, label_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Construct instance of huxley continuation problem from initial data.
prob = coco_prob();

% Set NTST size
prob = coco_set(prob, 'coll', 'NTST', 10);

% % Set step sizes
% h = 1e-1;
% prob = coco_set(prob, 'cont', 'h_min', h, 'h0', h, 'h_max', h);

% Set Continuation steps
PtMX = 250;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

prob = coco_set(prob, 'coll', 'MXCL', false);

% Set inital solution for trajectory segment
prob = ode_isol2coll(prob, 'homoclinic', func_list{:}, ...
                     single_segment.t_isol, single_segment.x_isol, single_segment.p0);

% Set initial solution of equilibrium point
prob = ode_isol2ep(prob, 'x_neg', func_list{:}, single_segment.x0_neg, ...
                   single_segment.p0);

% Glue that shit together, haumi ;)
prob = glue_single_segment(prob, data_lins, epsilon, eigvecs, eigvals);   

% Run COCO
coco(prob, run_new, [], 1, {'A', 'gamma', 'T', 'theta'}, p_range);
% coco(prob, run_new, [], 1, {'A', 'gamma', 'eps1', 'eps2', 'theta'}, p_range);

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
% Plot homoclinic method comparison diagram
compare_homoclinic_bifurcations(run_names, save_figure);

%-------------------%
%     Read Data     %
%-------------------%
% % Choose solution label
% % labels = coco_bd_labs(bdtest, 'ALL')
% label_plot = 1;
% 
% % Read solution
% [sol, ~] = coll_read_solution('homoclinic', run_new, label_plot);
% x_plot = sol.xbp;
% t_plot = sol.tbp;

% %--------------%
% %     Plot     %
% %--------------%
% % Create figure
% fig = figure(69); fig.Name = 'Test Homoclinic Continuation Plot'; clf;
% fig.Units = 'inches'; fig.Position = [3, 3, 12, 8]; fig.PaperSize = [12, 8];
% ax = gca();
% 
% hold(ax, 'on');
% 
% plot3(ax, x_plot(:, 1), x_plot(:, 2), x_plot(:, 3));
% 
% hold(ax, 'off')
% 
% % Axis limits
% % ax.XAxis.Limits = [-0.4, 0.3];
% % ax.YAxis.Limits = [-0.6, 0.25];
% % ax.ZAxis.Limits = [-0.2, 0.75];
% 
% % Axis Labels
% ax.XAxis.Label.String = '$G(t)$';
% ax.YAxis.Label.String = '$Q(t)$';
% ax.ZAxis.Label.String = '$I(t)$';
% 
% % Axis title
% title_str = sprintf('COCO Solution (run: $\\verb!%s!$, label: %d)', run_new, label_plot);
% ax.Title.String = title_str;
% 
% % Tick params
% ax.XAxis.TickDirection = 'in';
% ax.YAxis.TickDirection = 'in';
% ax.ZAxis.TickDirection = 'in';
% 
% % Figure stuff
% box(ax, 'on');
% grid(ax, 'on');
% 
% ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;
% view(45, 15.0);

function data_single_segment_out = merge_segments(run_in, label_in)
  % data_single_segment_out = merge_segments(run_in, label_in)
  %
  % Merge the unstable and stable segments into a single segment

  %-------------------%
  %     Read Data     %
  %-------------------%
  [sol1, ~]    = coll_read_solution('unstable', run_in, label_in);
  [sol2, ~]    = coll_read_solution('stable', run_in, label_in);
  [sol_neg, ~] = ep_read_solution('x0_neg', run_in, label_in);

  single_segment.p0 = sol_neg.p;
  single_segment.x0_neg = sol_neg.x;

  %--------------------%
  %     Merge Data     %
  %--------------------%
  % Grab state date and time data
  x1 = sol1.xbp(1:end-1, :);
  x2 = sol2.xbp(2:end, :);
  t1 = sol1.tbp(1:end-1);
  t2 = sol2.tbp(2:end);

  % Split data up
  G = [x1(:, 1); x2(2:end, 1)];
  Q = [x1(:, 2); x2(2:end, 2)];
  I = [x1(:, 3); x2(2:end, 3)];

  % Merge data
  single_segment.x_isol = [G, Q, I];
  single_segment.t_isol = [t1; t1(end) + t2(2:end)];

  %----------------%
  %     Output     %
  %----------------%
  data_single_segment_out = single_segment;

end

function prob_out = glue_single_segment(prob_in, data_in, epsilon_in, eigvecs_in, eigvals_in)

  % Set the COCO problem
  prob = prob_in;

  %-------------------------------%
  %     Read the Segment Data     %
  %-------------------------------%
  % Read segment data
  [data_coll, uidx_coll] = coco_get_func_data(prob, 'homoclinic.coll', 'data', 'uidx');
  [data_ep, uidx_ep] = coco_get_func_data(prob, 'x0.ep', 'data', 'uidx');
  % Grab the indices
  maps_coll = data_coll.coll_seg.maps;
  maps_ep = data_ep.ep_eqn;

  %-----------------------------------------------------%
  %     Glue Trajectory Segment Parameters Together     %
  %-----------------------------------------------------%
  % Glue parameters together
  prob = coco_add_glue(prob, 'parameters_shared', ...
                       uidx_coll(maps_coll.p_idx), uidx_ep(maps_ep.p_idx));

  %---------------------------------------------------------%
  %     Eigenvalue and Eigenvector Conditions: Unstable     %
  %---------------------------------------------------------%
  % Unstable eigenvector and eigenvalues
  vu = eigvecs_in{1};
  lu = eigvals_in{1};

  % Append boundary conditions to continue the eigenvalue and eigenvectors
  % of the Jacobian.
  % Unstable eigenvector
  prob = coco_add_func(prob, 'bcs_eig_unstable', @boundary_conditions_eig, [], ...
                       'zero', 'uidx', ...
                       [uidx_ep(maps_ep.x_idx); ...
                        uidx_coll(maps_coll.p_idx)], ...
                       'u0', [vu; lu]);

  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "vec_floquet" and "lam_floquet".
  uidx_eigu = coco_get_func_data(prob, 'bcs_eig_unstable', 'uidx');

  % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet
  vu_idx = [numel(uidx_eigu) - 3; numel(uidx_eigu) - 2; numel(uidx_eigu) - 1];
  lu_idx = numel(uidx_eigu);

  % Define active parameters for unstable eigenvector and eigenvalue
  prob = coco_add_pars(prob, 'eig_unstable', ...
                       [uidx_eigu(vu_idx); uidx_eigu(lu_idx)], ...
                       {'vu_1', 'vu_2', 'vu_3', 'lu'}, ...
                       'active');

  %---------------------------------------------------------%
  %     Eigenvalue and Eigenvector Conditions: Stable 1     %
  %---------------------------------------------------------%
  % Unstable eigenvector and eigenvalues
  vs1 = eigvecs_in{2};
  ls1 = eigvals_in{2};

  % Append boundary conditions to continue the eigenvalue and eigenvectors
  % of the Jacobian.
  % Unstable eigenvector
  prob = coco_add_func(prob, 'bcs_eig_stable1', @boundary_conditions_eig, [], ...
                       'zero', 'uidx', ...
                       [uidx_ep(maps_ep.x_idx); ...
                        uidx_coll(maps_coll.p_idx)], ...
                       'u0', [vs1; ls1]);

  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "vec_floquet" and "lam_floquet".
  uidx_eigs1 = coco_get_func_data(prob, 'bcs_eig_stable1', 'uidx');

  % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet
  vs1_idx = [numel(uidx_eigs1) - 3; numel(uidx_eigs1) - 2; numel(uidx_eigs1) - 1];
  ls1_idx = numel(uidx_eigs1);

  % Define active parameters for stable eigenvector 1 and eigenvalue 1
  prob = coco_add_pars(prob, 'eig_stable1', ...
                       [uidx_eigs1(vs1_idx); uidx_eigs1(ls1_idx)], ...
                       {'vs1_1', 'vs1_2', 'vs1_3', 'ls1'}, ...
                       'active');

  %---------------------------------------------------------%
  %     Eigenvalue and Eigenvector Conditions: Stable 2     %
  %---------------------------------------------------------%
  % Unstable eigenvector and eigenvalues
  vs2 = eigvecs_in{3};
  ls2 = eigvals_in{3};

  % Append boundary conditions to continue the eigenvalue and eigenvectors
  % of the Jacobian.
  % Unstable eigenvector
  prob = coco_add_func(prob, 'bcs_eig_stable2', @boundary_conditions_eig, [], ...
                      'zero', 'uidx', ...
                      [uidx_ep(maps_ep.x_idx); ...
                       uidx_coll(maps_coll.p_idx)], ...
                      'u0', [vs2; ls2]);

  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "vec_floquet" and "lam_floquet".
  uidx_eigs2 = coco_get_func_data(prob, 'bcs_eig_stable2', 'uidx');

  % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet
  vs2_idx = [numel(uidx_eigs2) - 3; numel(uidx_eigs2) - 2; numel(uidx_eigs2) - 1];
  ls2_idx = numel(uidx_eigs2);

  % Define active parameters for stable eigenvector 2 and eigenvalue 2
  prob = coco_add_pars(prob, 'eig_stable2', ...
                      [uidx_eigs2(vs2_idx); uidx_eigs2(ls2_idx)], ...
                      {'vs2_1', 'vs2_2', 'vs2_3', 'ls'}, ...
                      'active');

  %-------------------------------------%
  %     Initial Boundary Conditions     %
  %-------------------------------------%
  % Apply boundary conditions for the end points
  prob = coco_add_func(prob, 'bcs_homo', @boundary_conditions_initial, [], 'zero', 'uidx', ...
                      [uidx_coll(maps_coll.x0_idx); ...
                       uidx_coll(maps_coll.x1_idx); ...
                       uidx_ep(maps_ep.x_idx); ...
                       uidx_eigu(vu_idx);
                       uidx_eigs1(vs1_idx);
                       uidx_eigs2(vs2_idx)], ...
                      'u0', epsilon_in);

  % Get u-vector indices from function
  uidx_eps = coco_get_func_data(prob, 'bcs_homo', 'uidx');
  % Get epsilon indices
  epsilon_idx = [numel(uidx_eps) - 2; numel(uidx_eps) - 1; numel(uidx_eps)];

  %-----------------------------------%
  %     Define Problem Parameters     %
  %-----------------------------------%
  % Define system parameters
  prob = coco_add_pars(prob, 'pars_system', ...
                      uidx_coll(maps_coll.p_idx), ...
                      {data_in.pnames{:}}, ...
                      'inactive');

  % Define problem parameters
  prob = coco_add_pars(prob, 'pars_eps', ...
                       uidx_eps(epsilon_idx), {'eps1', 'eps2', 'theta'}, ...
                       'inactive');

  % Define trajectory periods
  prob = coco_add_pars(prob, 'pars_T', ...
                       uidx_coll(maps_coll.T_idx), 'T', ...
                       'inactive');

  %----------------%
  %     Output     %
  %----------------%
  % Output problem structure
  prob_out = prob;

end