%-------------------------------------------------------------------------%
%%                        Hopf to Periodic Orbit                         %%
%-------------------------------------------------------------------------%
% Continue a family of periodic orbits emanating from the Hopf
% bifurcation with 'ode_HB2po'.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.hopf_to_PO;
% Which run this continuation continues from
run_old = run_names.move_hopf;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'H_PT');

% Print to console
fprintf("~~~ Initial Periodic Orbit: Fifth Run (hopf_to_PO.m) ~~~ \n");
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
prob = coco_set(prob, 'coll', 'NTST', 60);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 800;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% % Set step sizes
% h_size = 1e0;
% prob = coco_set(prob, 'cont', 'h_min', h_size, 'h0', h_size, 'h_max', h_size);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Turn off bifurcation detection
% prob = coco_set(prob, 'po', 'bifus', 'off');

% Continue from Hopf bifurcation
prob = ode_HB2po(prob, '', run_old, label_old);
% prob = ode_HB2po(prob, '', run_old, label_old, '-var', eye(3));

% Follow non trivial solutions
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   xpos, sol.p);
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, ...
                   xneg, sol.p);
prob = ode_isol2ep(prob, 'x0',   funcs.field{:}, ...
                   x0,   sol.p);

% Glue parameters
prob = glue_parameters(prob);

% Event for Nonlinear Photonics abstract
prob = coco_add_event(prob, 'PO_PT', 'gamma', 3.54e-2);

% Halfway down from \gamma = 0 to the H_PT
% prob = coco_add_event(prob, 'PO_PT', 'gamma', 6.1953e-2);

% Run COCO
bd_PO = coco(prob, run_new, [], 1, {'gamma', 'A'}, gamma_range);

%-------------------------------------------------------------------------%
%%                            Testing Things                             %%
%-------------------------------------------------------------------------%
% Solution to plot
label_plot = sort(coco_bd_labs(coco_bd_read(run_new), 'PO_PT'));
label_plot = label_plot(1);

% label_plot = sort(coco_bd_labs(coco_bd_read(run_new), 'EP'));
% label_plot = max(label_plot);

% Create plots
plot_hopf_to_PO_solution(run_new, label_plot);

%------------------------------------------------%
%     Calculate Eigenvalues and Eigenvectors     %
%------------------------------------------------%
% % Read one of the solutions
% chart = coco_read_solution('po.orb.coll.var', run_new, label_plot, 'chart');
% data  = coco_read_solution('po.orb.coll', run_new, label_plot, 'data');
% 
% % Create monodrony matrix
% M1 = chart.x(data.coll_var.v1_idx);
% 
% fprintf('~~~ Monodromy Matrix ~~~\n');
% fprintf('(%.7f, %.7f, %.7f)\n', M1(1, :));
% fprintf('(%.7f, %.7f, %.7f)\n', M1(2, :));
% fprintf('(%.7f, %.7f, %.7f)\n\n', M1(3, :));
% 
% % Calculate eigenvalues and eigenvectors
% [v, d] = eig(M1);
% 
% % Find index for stable eigenvector? < 1
% ind = find(abs(diag(d)) < 1);
% 
% % Stable eigenvector
% vec0 = -v(:, ind);
% % Stable eigenvalue (Floquet thingie)
% lam0 = d(ind, ind);
% 
% % Do all the same but with the function I defined
% % [vec0, lam0] = calculate_stable_floquet(run_new, label_plot);
% 
% fprintf('\n~~~ Eigenvector and Eigenvalue ~~~\n');
% fprintf('vec0 (numeric)  = (%f, %f, %f) \n', vec0);
% fprintf('lam0 (numeric)  = %s \n\n', lam0);


%-------------------------------------------------------------------------%
%%                               FUNCTIONS                               %%
%-------------------------------------------------------------------------%
function prob_out = glue_parameters(prob_in)
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
  % Read index data periodic orbit segment
  [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');

  % Read index data for the variational problem
  % [data_var, uidx_var] = coco_get_func_data(prob, 'po.orb.coll.var', 'data', 'uidx');

  % Read index data equilibrium points
  [data1, uidx1] = coco_get_func_data(prob, 'xpos.ep', 'data', 'uidx');
  [data2, uidx2] = coco_get_func_data(prob, 'xneg.ep', 'data', 'uidx');
  [data3, uidx3] = coco_get_func_data(prob, 'x0.ep',   'data', 'uidx');

  % Index mapping
  maps = data.coll_seg.maps;
  % maps_var = data_var.coll_var;
  maps1 = data1.ep_eqn;
  maps2 = data2.ep_eqn;
  maps3 = data3.ep_eqn;

  %-------------------------%
  %     Glue Parameters     %
  %-------------------------%
  prob = coco_add_glue(prob, 'glue_p1', uidx(maps.p_idx), uidx1(maps1.p_idx));
  prob = coco_add_glue(prob, 'glue_p2', uidx(maps.p_idx), uidx2(maps2.p_idx));
  prob = coco_add_glue(prob, 'glue_p3', uidx(maps.p_idx), uidx3(maps3.p_idx));

  %------------------------%
  %     Add Parameters     %
  %------------------------%
  % % Add parameters for each component of the monodromy matrix
  % prob = coco_add_pars(prob, 'pars_var', ...
  %                      uidx_var(maps_var.v0_idx,:), ...
  %                      {'s1', 's2', 's3', ...
  %                       's4', 's5', 's6', ...
  %                       's7', 's8', 's9'});

  %----------------------------------------%
  %     Eigenvalue Boundary Conditions     %
  %----------------------------------------%
  % % % Data structure containing the toolbox id for the monodromy matrix data
  % % data_eig.tbid = 'hopf_po.po.orb.coll';
  % 
  % % Apply the boundary conditions for the eigenvalues and eigenvectors of the
  % % monodromy matrix
  % prob = coco_add_func(prob, 'bcs_eigen', @boundary_conditions_eig, data_var, ...
  %                      'zero', 'uidx', ...
  %                      [uidx_var(maps_var.v1_idx(:, 1)); ...
  %                       uidx_var(maps_var.v1_idx(:, 2)); ...
  %                       uidx_var(maps_var.v1_idx(:, 3))], ...
  %                      'u0', [vec_floquet_in; lam_floquet_in]);
  % 
  % % Get u-vector indices from this coco_add_func call, including the extra
  % % indices from the added "vec_floquet" and "lam_floquet".
  % uidx_eig = coco_get_func_data(prob, 'bcs_eigen', 'uidx');
  % 
  % % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet]
  % data_out.vec_floquet_idx = uidx_eig(end-3:end-1);
  % data_out.lam_floquet_idx = uidx_eig(end);

  %----------------%
  %     Output     %
  %----------------%
  prob_out = prob;

end

function [data_in, y_out] = boundary_conditions_eig(prob_in, data_in, u_in)
  % [data_in, y_out] = boundary_conditions_eig(prob_in, data_in, u_in)
  % 
  % COCO compatible encoding for the boundary conditions of the eigenvalues and
  % eigenvectors of the monodromy matrix. Ensures they are eigenvectors and
  % values of the monodromy matrix, and ensures that the eigenvector is
  % normalised.
  %
  % Input
  % ----------
  % prob_in : COCO problem structure
  %     Continuation problem structure.
  % data_in : structure
  %     Problem data structure contain with function data.
  % u_in : array (floats?)
  %     Total u-vector of the continuation problem. This function
  %     only utilises the following (as imposed by coco_add_func):
  %          * u_in(1:9)   - The monodromy matrix,
  %          * u_in(10:12) - The stable Floquet vector (vec_floquet),
  %          * u_in(13)    - The stable Floquet multiplier (lam_floquet).
  %
  % Output
  % ----------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Not actually output here but you need to have it for COCO.

  xdim = data_in.xdim;

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Monodromy matrix indices
  v1_idx      = u_in(1 : (xdim ^ 2));
  % Eigenvector
  vec_floquet = u_in(end-xdim : end-1);
  % Eigenvalue
  lam_floquet = u_in(end);

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Reshape indices for monodromy matrix 
  M1 = reshape(v1_idx, xdim, xdim);

  % fprintf('~~~ Monodromy Matrix ~~~\n');
  % fprintf('(%.7f, %.7f, %.7f)\n', M1(1, :));
  % fprintf('(%.7f, %.7f, %.7f)\n', M1(2, :));
  % fprintf('(%.7f, %.7f, %.7f)\n\n', M1(3, :));

  %---------------------------------------%
  %     Calculate Boundary Conditions     %
  %---------------------------------------%
  % Eigenvalue equations
  eig_eqn = (M1 * vec_floquet) - (lam_floquet * vec_floquet);

  % Unit vector equations
  vec_eqn = (vec_floquet' * vec_floquet) - 1;

  %----------------%
  %     Output     %
  %----------------%
  % Boundary conditions
  y_out = [eig_eqn;
           vec_eqn];

end