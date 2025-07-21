function [data_in, y_out] = bcs_Wsq_initial(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_Wsq_initial(prob_in, data_in, u_in)
  %
  % COCO compatible encoding for the "final" boundary conditions of the
  % trajectory segment of the stable manifold of x_pos. The segment
  % will intersect a defined plane, \Sigma, 
  %           x_Wsq(0) \in \Sigma ,
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
  %          * u_in(1:3) - The initial point of manifold trajectory segment (x0_Wsq),
  %          * u_in(4:6) - The initial point of \gamma_{\vartheta}} (x0_seg).
  %
  % Output
  % ----------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Not actually output here but you need to have it for COCO.

  % State space dimension
  xdim = data_in.xdim;
  
  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Initial vector of segment 1
  x0_Wsq = u_in(1:3);
  % Initial vector of segment 2
  x0_seg = u_in(4:6);

  %----------------%
  %     Output     %
  %----------------%
  % Boundary conditions
  y_out = x0_Wsq(2) - x0_seg(2);

end
