function [data_in, y_out] = monitor_q_distance(prob_in, data_in, u_in)
  % [data_in, y_out] = monitor_q_distance(prob_in, data_in, u_in)
  % 
  % Monitor the distance in the (X, Z) or (G, I) plane from the point
  % \gamma_{\vartheta_{old}} to the equilibrium point q.
  %
  % Parameters
  % ----------
  % prob_in : COCO problem structure
  %     Continuation problem structure.
  % data_in : structure
  %     Problem data structure contain with function data.
  % u_in : array (floats?)
  %     Total u-vector of the continuation problem. This function
  %     only utilises the following (as imposed by coco_add_func):
  %          * u_in(1:3) - Initial point of \gamma_{\vartheta}
  %          * u_in(4:6) - Equilibrium point q (x_{+}}.
  %
  % Returns
  % -------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Function data structure to give dimensions of parameter and state
  %     space.

  % Original vector field dimensions
  xdim = data_in.xdim;

  %---------------%
  %     Input     %
  %---------------%
  % Initial point of the periodic orbit
  x0_seg = u_in(1 : xdim);
  % Final point of the periodic orbit
  xpos   = u_in(xdim+1 : 2*xdim);

  %----------------%
  %     Output     %
  %----------------%
  y_out = x0_seg(2) - xpos(2);

end