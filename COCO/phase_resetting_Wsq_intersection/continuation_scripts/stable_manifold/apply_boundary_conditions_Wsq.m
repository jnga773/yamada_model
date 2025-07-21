function prob_out = apply_boundary_conditions_Wsq(prob_in, data_in, bcs_funcs_in, epsilon_in)
  % prob_out = apply_boundary_conditions_Wsq(prob_in, data_in, bcs_funcs_in, epsilon_in)
  %
  % This function reads index data for the stable periodic orbit segment and equilibrium points,
  % glues the COLL and EP parameters together, applies periodic orbit boundary conditions,
  % and adds variational problem matrix parameters.
  %
  % Parameters
  % ----------
  % prob_in : COCO problem structure
  %     Input continuation problem structure.
  % bcs_funcs_in : List of functions
  %     Structure containing boundary condition functions.
  % data_in : structure
  %     Problem data structure containing initialisation information
  %     for the segments.
  % epsilon_in : double
  %     Epsilon parameter for boundary conditions.
  %
  % Returns
  % -------
  % prob_out : COCO problem structure
  %     Output continuation problem structure with applied boundary conditions.
  %
  % See Also
  % --------
  % coco_get_func_data, coco_add_glue, coco_add_func, coco_add_pars

  %---------------%
  %     Input     %
  %---------------%
  % Set the COCO problem
  prob = prob_in;

  % Original vector field dimensions
  xdim            = data_in.xdim;
  pdim            = data_in.pdim;
  % Create data structure for original vector field dimensions
  dim_data.xdim   = xdim;
  dim_data.pdim   = pdim;
  dim_data.p_maps = data_in.p_maps;

  %-------------------------------%
  %     Read the Segment Data     %
  %-------------------------------%
  % Extract toolbox data and context-dependent index set for each of the orbit
  % segments.
  [data1, uidx1]   = coco_get_func_data(prob, 'seg1.coll', 'data', 'uidx');
  [data2, uidx2]   = coco_get_func_data(prob, 'seg2.coll', 'data', 'uidx');

  % Grab the indices from each of the orbit segments
  maps1 = data1.coll_seg.maps;
  maps2 = data2.coll_seg.maps;

  % Read index data for equilibrium points
  [data_pos, uidx_pos] = coco_get_func_data(prob, 'xpos.ep', 'data', 'uidx');
  % Index mapping
  maps_pos   = data_pos.ep_eqn;

  % Extract toolbox data and context-dependent index set for manifold segments
  [dataW, uidxW] = coco_get_func_data(prob, 'Wsq.coll', 'data', 'uidx');
  % Index mapping
  mapsW = dataW.coll_seg.maps;

  %-----------------------------------------------------%
  %     Glue Trajectory Segment Parameters Together     %
  %-----------------------------------------------------%
  % "Glue" segment periods together
  prob = coco_add_glue(prob, 'glue_T', uidx1(maps1.T_idx), uidx2(maps2.T_idx));
                       
  % "Glue" segment parameters together
  prob = coco_add_glue(prob, 'glue_pars', uidx1(maps1.p_idx), uidx2(maps2.p_idx));

  % Glue equilibrium point and segment parameters together
  prob = coco_add_glue(prob, 'glue_ep', ...
                       uidx1(maps1.p_idx(1:pdim)), uidx_pos(maps_pos.p_idx));

  % Glue phase resetting and manifold segment parameters together
  prob = coco_add_glue(prob, 'glue_Wsq', ...
                       uidx1(maps1.p_idx(1:pdim)), uidxW(mapsW.p_idx));

  %---------------------------------------------%
  %     Phase Resetting Boundary Conditions     %
  %---------------------------------------------%
  % Boundary condition function list
  bcs_PR = bcs_funcs_in.bcs_PR;

  % Add boundary conditions for four segments
  prob = coco_add_func(prob, 'bcs_PR', bcs_PR{:}, dim_data, 'zero', 'uidx', ...
                       [uidx1(maps1.x0_idx);
                        uidx2(maps2.x0_idx);
                        uidx1(maps1.x1_idx);
                        uidx2(maps2.x1_idx);
                        uidx1(maps1.p_idx)]);

  %-------------------------------------%
  %     Initial Boundary Conditions     %
  %-------------------------------------%
  % Boundary condition function list
  bcs_Wsq_initial = bcs_funcs_in.bcs_Wsq_initial;

  % Append hyperplane conditions with parameters 'seg_u' and 'seg_s' for the unstable
  % and stable segments, respectively.
  prob = coco_add_func(prob, 'bcs_Wsq_initial', bcs_Wsq_initial{:}, data1, ...
                       'inactive', 'Wsq_dist', 'uidx', ...
                       [uidxW(mapsW.x0_idx); ...
                        uidx2(maps2.x0_idx)]);

  %-----------------------------------%
  %     Final Boundary Conditions     %
  %-----------------------------------%
  % Boundary condition function list
  bcs_Wsq_final = bcs_funcs_in.bcs_Wsq_final;

  % Apply the boundary conditions for the initial points near the equilibrium.
  % Here, epsilon_in is appended as an input to the u-vector input for the function
  % @bcs_initial.
  prob = coco_add_func(prob, 'bcs_Wsq_final', bcs_Wsq_final{:}, data_pos.data, 'zero', 'uidx', ...
                       [uidxW(mapsW.x1_idx); ...
                        uidx_pos(maps_pos.x_idx)], ...
                       'u0', epsilon_in);

  % Get u-vector indices from this coco_add_func call, including the extra
  % indices from the added "eps"
  uidx_eps = coco_get_func_data(prob, 'bcs_Wsq_final', 'uidx');
  
  % Grab eigenvector and value indices from u-vector [eps]
  data_out.eps_idx = numel(uidx_eps);
  
  % Save data
  prob = coco_add_slot(prob, 'bcs_Wsq_final', @coco_save_data, data_out, 'save_full');

  %-----------------------------------%
  %     Define Problem Parameters     %
  %-----------------------------------%
  % Monitor distance from periodic orbit to q in (G, I) plane
  prob = coco_add_func(prob, 'monitor_q', @monitor_q_distance, dim_data, ...
                       'active', 'q_dist', 'uidx', ...
                       [uidx2(maps2.x0_idx); uidx_pos(maps_pos.x_idx)]);
                       
  % Define trajectory periods
  prob = coco_add_pars(prob, 'pars_T', ...
                       uidxW(mapsW.T_idx), 'T_Wsq', 'inactive');

  % Define epsilon parameters
  prob = coco_add_pars(prob, 'pars_eps', ...
                       uidx_eps(data_out.eps_idx), 'epsilon', 'inactive');

  % Monitor distance from periodic orbit and (G, I) plane of manifold
  % segment.
  prob = coco_add_func(prob, 'monitor_q_theta_old', @monitor_q_perturb, dim_data, ...
                       'active', {'A_perturb', 'theta_perturb'}, 'uidx', ...
                       [uidx2(maps2.x0_idx); uidxW(mapsW.x0_idx)]);

  %----------------%
  %     Output     %
  %----------------%
  % Output problem structure
  prob_out = prob;

end