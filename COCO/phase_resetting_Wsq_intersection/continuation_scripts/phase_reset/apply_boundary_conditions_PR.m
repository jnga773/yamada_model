function prob_out = apply_boundary_conditions_PR(prob_in, data_in, bcs_funcs_in)
  % prob_out = apply_boundary_conditions_PR(prob_in)
  %
  % Applies the various boundary conditions, adds monitor functions
  % for the singularity point, glue parameters and other things
  % together, and also adds some user defined points. Also
  % applied some COCO settings to the continuation problem.
  %
  % Parameters
  % ----------
  % prob_in : COCO problem structure
  %     Continuation problem structure.
  % data_in : structure
  %     Problem data structure containing initialisation information
  %     for the segments.
  % bcs_funcs_in : list of functions
  %     List of all of the boundary condition functions for each
  %     phase resetting segment.
  %
  % Returns
  % -------
  % prob_out : COCO problem structure
  %     Continuation problem structure.

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

  % Extract toolbox data and mappings for the equilibrium point
  [data_pos, uidx_pos] = coco_get_func_data(prob, 'xpos.ep', 'data', 'uidx');
  maps_pos = data_pos.ep_eqn;

  %----------------------------------------%
  %     Glue Trajectory Segment Things     %
  %----------------------------------------%
  % "Glue" segment periods together
  prob = coco_add_glue(prob, 'glue_T', uidx1(maps1.T_idx), uidx2(maps2.T_idx));
                       
  % "Glue" segment parameters together
  prob = coco_add_glue(prob, 'glue_pars', uidx1(maps1.p_idx), uidx2(maps2.p_idx));

  % Glue equilibrium point and segment parameters together
  prob = coco_add_glue(prob, 'glue_ep', uidx1(maps1.p_idx(1:pdim)), uidx_pos(maps_pos.p_idx));

  %---------------------------------%
  %     Add Boundary Conditions     %
  %---------------------------------%
  % Boundary condition function list
  bcs_PR = bcs_funcs_in.bcs_PR;

  % Add boundary conditions for four segments
  prob = coco_add_func(prob, 'bcs_PR', bcs_PR{:}, dim_data, 'zero', 'uidx', ...
                       [uidx1(maps1.x0_idx);
                        uidx2(maps2.x0_idx);
                        uidx1(maps1.x1_idx);
                        uidx2(maps2.x1_idx);
                        uidx1(maps1.p_idx)]);
                        
  %------------------------%
  %     Add Parameters     %
  %------------------------%
  % Add parameter to monitor I at \gamma_{\theta_{o}}
  prob = coco_add_pars(prob, 'pars_I_theta_n', ...
                       uidx2(maps2.x0_idx(3)), 'I_theta_n', ...
                       'active');

  % Monitor distance from periodic orbit to q in (G, I) plane
  prob = coco_add_func(prob, 'monitor_q', @monitor_q_distance, dim_data, ...
                       'active', 'q_dist', 'uidx', ...
                       [uidx2(maps2.x0_idx); uidx_pos(maps_pos.x_idx)]);

  %----------------%
  %     Output     %
  %----------------%
  % Output problem structure
  prob_out = prob;

end
