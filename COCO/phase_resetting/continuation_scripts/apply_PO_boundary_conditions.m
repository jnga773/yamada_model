function prob_out = apply_PO_boundary_conditions(prob_in, bcs_PO_in)
  % prob_out = apply_PO_boundary_conditions(prob_in)
  %
  % Glue the COLL and EP parameters together, and apply the 'new'
  % periodic orbit boundary conditions.

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
  maps  = data.coll_seg.maps;
  maps1 = data1.ep_eqn;
  maps2 = data2.ep_eqn;
  maps3 = data3.ep_eqn;

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