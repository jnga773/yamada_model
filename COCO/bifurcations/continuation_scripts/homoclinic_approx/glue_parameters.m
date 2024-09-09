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
  [data_PO, uidx_PO] = coco_get_func_data(prob, 'homo.po.orb.coll', 'data', 'uidx');

  % Read index data equilibrium points
  [data_x0, uidx_x0] = coco_get_func_data(prob, 'x0.ep', 'data', 'uidx');

  % Index mapping
  maps_PO = data_PO.coll_seg.maps;
  maps_x0 = data_x0.ep_eqn;

  %-------------------------%
  %     Glue Parameters     %
  %-------------------------%
  prob = coco_add_glue(prob, 'shared_parameters', ...
                       uidx_PO(maps_PO.p_idx), uidx_x0(maps_x0.p_idx));

  %----------------%
  %     Output     %
  %----------------%
  prob_out = prob;

end