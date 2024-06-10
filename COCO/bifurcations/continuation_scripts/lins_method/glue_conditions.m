function prob_out = glue_conditions(prob_in, data_in, epsilon_in, eigvecs_in, eigvals_in)
  % prob_out = glue_conditions(prob_in, data_in, epsilon_in)
  %
  % Encoding of the initial and final boundary conditions of the two trajectory segments.
  % 
  % The boundary conditions are applied as added monitor functions,
  % with inactive parameters 'seg_u' and 'seg_s', which monitor the
  % distance of the trajectory segment from the defined point on the
  % \Sigma plane, for the unstable and stable manifolds, respectively.
  %
  % The system parameters are also added, with the coco_glue function
  % setting the parameters of each segment to be the same.
  % 
  % The epsilon distance parameters are then saved into the data_in
  % structure such that they can be called in sucessive runs.
  %
  % Finally, COCO events are added for when the trajectory segments
  % hit the \Sigma plane, i.e., when 'seg_u' or 'seg_s' = 0.
  %
  % Input
  % ----------
  % prob_in : COCO problem structure
  %     Continuation problem structure.
  % data_in : structure
  %     Problem data structure containing boundary condition information.
  % epsilon_in : array (floats)
  %     Array of separations at end points from orbit and equilibrium.
  %
  % Output
  % ----------
  % prob_out        : COCO problem structure
  %     Continuation problem structure.

  % Set the COCO problem
  prob = prob_in;

  %-------------------------------%
  %     Read the Segment Data     %
  %-------------------------------%
  % Extract toolbox data and context-dependent index set for the unstable
  % and stable manfiolds.
  [data_u, uidx_u] = coco_get_func_data(prob, 'unstable.coll', 'data', 'uidx');
  [data_s, uidx_s] = coco_get_func_data(prob, 'stable.coll', 'data', 'uidx');

  % Grab the indices from data_u and data_s
  maps_u = data_u.coll_seg.maps;
  maps_s = data_s.coll_seg.maps;

  % uidx_u(maps_u.x0_idx) - u-vector indices for the initial point of the
  %                         unstable manifold trajectory segment - x_u(0).
  % uidx_u(maps_u.x1_idx) - u-vector indices for the final point of the
  %                         unstable manifold trajectory segment - x_u(T1).

  % uidx_s(maps_s.x0_idx) - u-vector indices for the initial point of the
  %                         stable manifold trajectory segment - x_s(0).
  % uidx_s(maps_s.x1_idx) - u-vector indices for the final point of the
  %                         stable manifold trajectory segment - x_s(T2).

  % uidx_u(maps_u.p_idx)  - u-vector indices for the system parameters.
  % uidx_s(maps_s.p_idx)  - u-vector indices for the system parameters.
  % uidx_u(maps_u.T_idx)  - u-vector index for the period of the unstable
  %                         manifold trajectory segment.
  % uidx_s(maps_s.T_idx)  - u-vector index for the period of the stable
  %                         manifold trajectory segment.

  % Extract toolbox data and indices for the x_neg and x_pos equilibrium points
  % from the 'ep' toolbox structures.
  [data_neg, uidx_neg] = coco_get_func_data(prob, 'x_neg.ep', 'data', 'uidx');
  [data_pos, uidx_pos] = coco_get_func_data(prob, 'x_pos.ep', 'data', 'uidx');

  % Grab the indices from data_neg and data_pos
  maps_neg = data_neg.ep_eqn;
  maps_pos = data_pos.ep_eqn;

  % uidx_neg(maps_neg.x_idx) - u-vector indices for the x_neg equilibrium point.
  % uidx_pos(maps_pos.x_idx) - u-vector indices for the x_pos equilibrium point.

  % uidx_neg(maps_neg.p_idx) - u0-vector indices for the system parameters
  % uidx_pos(maps_pos.p_idx) - u0-vector indices for the system parameters

  %-----------------------------------------------------%
  %     Glue Trajectory Segment Parameters Together     %
  %-----------------------------------------------------%
  % All segments have the same system parameters, so "glue" them together,
  % i.e., let COCO know that they are the same thing.
  prob = coco_add_glue(prob, 'parameters_shared_coll', ...
                       uidx_u(maps_u.p_idx), uidx_s(maps_s.p_idx));
  prob = coco_add_glue(prob, 'parameters_shared_x_neg', ...
                       uidx_u(maps_u.p_idx), uidx_neg(maps_neg.p_idx));
  prob = coco_add_glue(prob, 'parameters_shared_x_pos', ...
                       uidx_u(maps_u.p_idx), uidx_pos(maps_pos.p_idx));

  %---------------------------------------------------------%
  %     Eigenvalue and Eigenvector Conditions: Unstable     %
  %---------------------------------------------------------%
  % Unstable eigenvector and eigenvalues
  vu = data_in.eigvecs{1};
  lu = data_in.eigvals{1};

  % Append boundary conditions to continue the eigenvalue and eigenvectors
  % of the Jacobian.
  % Unstable eigenvector
  prob = coco_add_func(prob, 'bcs_eig_unstable', @bcs_eig, data_u, ...
                       'zero', 'uidx', ...
                       [uidx_neg(maps_neg.x_idx); ...
                        uidx_u(maps_u.p_idx)], ...
                       'u0', [vu; lu]);
  
  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "vec_floquet" and "lam_floquet".
  uidx_eigu = coco_get_func_data(prob, 'bcs_eig_unstable', 'uidx');
  
  % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet
  % data_out.vu_idx = uidx_eigu(end-3:end-1);
  % data_out.lu_idx = uidx_eigu(end);
  data_out.vu_idx = [numel(uidx_eigu) - 3; numel(uidx_eigu) - 2; numel(uidx_eigu) - 1];
  data_out.lu_idx = numel(uidx_eigu);

  % Save data structure to be called later with coco_get_func_data 
  prob = coco_add_slot(prob, 'bcs_eig_unstable', @coco_save_data, data_out, 'save_full');

  % Define active parameters for unstable eigenvector and eigenvalue
  prob = coco_add_pars(prob, 'eig_unstable', ...
                       [uidx_eigu(data_out.vu_idx); uidx_eigu(data_out.lu_idx)], ...
                       {'vu_1', 'vu_2', 'vu_3', 'lu'}, ...
                       'active');

  %---------------------------------------------------------%
  %     Eigenvalue and Eigenvector Conditions: Stable 1     %
  %---------------------------------------------------------%
  % Unstable eigenvector and eigenvalues
  vs1 = data_in.eigvecs{2};
  ls1 = data_in.eigvals{2};

  % Append boundary conditions to continue the eigenvalue and eigenvectors
  % of the Jacobian.
  % Unstable eigenvector
  prob = coco_add_func(prob, 'bcs_eig_stable1', @bcs_eig, data_u, ...
                       'zero', 'uidx', ...
                       [uidx_neg(maps_neg.x_idx); ...
                        uidx_u(maps_u.p_idx)], ...
                       'u0', [vs1; ls1]);
  
  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "vec_floquet" and "lam_floquet".
  uidx_eigs1 = coco_get_func_data(prob, 'bcs_eig_stable1', 'uidx');
  
  % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet
  % data_out.vs1_idx = uidx_eigs1(end-3:end-1);
  % data_out.ls1_idx = uidx_eigs1(end);
  data_out.vs1_idx = [numel(uidx_eigs1) - 3; numel(uidx_eigs1) - 2; numel(uidx_eigs1) - 1];
  data_out.ls1_idx = numel(uidx_eigs1);

  % Save data structure to be called later with coco_get_func_data 
  prob = coco_add_slot(prob, 'bcs_eig_stable1', @coco_save_data, data_out, 'save_full');

  % Define active parameters for stable eigenvector 1 and eigenvalue 1
  prob = coco_add_pars(prob, 'eig_stable1', ...
                       [uidx_eigs1(data_out.vs1_idx); uidx_eigs1(data_out.ls1_idx)], ...
                       {'vs1_1', 'vs1_2', 'vs1_3', 'ls1'}, ...
                       'active');

  %---------------------------------------------------------%
  %     Eigenvalue and Eigenvector Conditions: Stable 2     %
  %---------------------------------------------------------%
  % Unstable eigenvector and eigenvalues
  vs2 = data_in.eigvecs{3};
  ls2 = data_in.eigvals{3};

  % Append boundary conditions to continue the eigenvalue and eigenvectors
  % of the Jacobian.
  % Unstable eigenvector
  prob = coco_add_func(prob, 'bcs_eig_stable2', @bcs_eig, data_u, ...
                       'zero', 'uidx', ...
                       [uidx_neg(maps_neg.x_idx); ...
                        uidx_u(maps_u.p_idx)], ...
                       'u0', [vs2; ls2]);
  
  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "vec_floquet" and "lam_floquet".
  uidx_eigs2 = coco_get_func_data(prob, 'bcs_eig_stable2', 'uidx');
  
  % Grab eigenvector and value indices from u-vector [vec_floquet, lam_floquet
  % data_out.vs2_idx = uidx_eigs2(end-3:end-1);
  % data_out.ls2_idx = uidx_eigs2(end);
  data_out.vs2_idx = [numel(uidx_eigs2) - 3; numel(uidx_eigs2) - 2; numel(uidx_eigs2) - 1];
  data_out.ls2_idx = numel(uidx_eigs2);

  % Save data structure to be called later with coco_get_func_data 
  prob = coco_add_slot(prob, 'bcs_eig_stable2', @coco_save_data, data_out, 'save_full');

  % Define active parameters for stable eigenvector 2 and eigenvalue 2
  prob = coco_add_pars(prob, 'eig_stable2', ...
                       [uidx_eigs2(data_out.vs2_idx); uidx_eigs2(data_out.ls2_idx)], ...
                       {'vs2_1', 'vs2_2', 'vs2_3', 'ls'}, ...
                       'active');

  %-------------------------------------%
  %     Initial Boundary Conditions     %
  %-------------------------------------%
  % Apply the boundary conditions for the initial points near the equilibrium.
  % Here, epsilon_in is appended as an input to the u-vector input for the function
  % @bcs_initial.
  prob = coco_add_func(prob, 'bcs_initial', @bcs_initial, data_u, 'zero', 'uidx', ...
                       [uidx_u(maps_u.x0_idx); ...
                        uidx_s(maps_s.x1_idx); ...
                        uidx_neg(maps_neg.x_idx); ...
                        uidx_eigu(data_out.vu_idx);
                        uidx_eigs1(data_out.vs1_idx);
                        uidx_eigs2(data_out.vs2_idx)], ...
                       'u0', data_in.epsilon);

  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "epsilon_in".
  uidx_eps = coco_get_func_data(prob, 'bcs_initial', 'uidx');

  % Grab epsilon parameters indices from u-vector [eps1, eps2, theta]
  % data_out.epsilon_idx = uidx_eps(end-2:end);
  data_out.epsilon_idx = [numel(uidx_eps) - 2; numel(uidx_eps) - 1; numel(uidx_eps)];

  % Also save pt0 and normal vector here
  data_out.pt0    = data_in.pt0;
  data_out.normal = data_in.normal;

  % Save data structure to be called later with coco_get_func_data 
  prob = coco_add_slot(prob, 'bcs_initial', @coco_save_data, data_out, 'save_full');

  %-----------------------------------%
  %     Final Boundary Conditions     %
  %-----------------------------------%
  % Append hyperplane conditions with parameters 'seg_u' and 'seg_s' for the unstable
  % and stable segments, respectively.
  data_in.xdim = data_u.xdim;
  prob = coco_add_func(prob, 'bcs_final', @bcs_final, data_in, ...
                       'inactive', {'seg_u', 'seg_s'}, 'uidx', ...
                       [uidx_u(maps_u.x1_idx); ...
                        uidx_s(maps_s.x0_idx); ...
                        uidx_pos(maps_pos.x_idx)]);

  %-----------------------------------%
  %     Define Problem Parameters     %
  %-----------------------------------%
  % Define epsilon parameters
  prob = coco_add_pars(prob, 'pars_eps', ...
                       uidx_eps(data_out.epsilon_idx), {'eps1', 'eps2', 'theta'}, ...
                       'inactive');

  % Define trajectory periods
  prob = coco_add_pars(prob, 'pars_T', ...
                       [uidx_u(maps_u.T_idx); uidx_s(maps_s.T_idx)], {'T1', 'T2'}, ...
                       'inactive');

  %--------------------------------------------%
  %     Add Event for Hitting \Sigma Plane     %
  %--------------------------------------------%
  % Add event for when trajectory hits \Sigma plane for unstable and stable manifolds.
  prob = coco_add_event(prob, 'DelU', 'seg_u', 0);
  prob = coco_add_event(prob, 'DelS', 'seg_s', 0);

  %----------------%
  %     Output     %
  %----------------%
  % Save data structure to be called later with coco_get_func_data 
  % prob = coco_add_slot(prob, 'apply_bcs', @coco_save_data, data_out, 'save_full');
  % prob = coco_add_slot(prob, 'bcs_initial', @coco_save_data, data_out, 'save_full');

  % Output problem structure
  prob_out = prob;

end