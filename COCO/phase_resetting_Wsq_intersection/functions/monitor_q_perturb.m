function [data_in, y_out] = monitor_q_perturb(prob_in, data_in, u_in)
  % [data_in, y_out] = monitor_q_distance(prob_in, data_in, u_in)
  % 
  % Monitor the distance in the (X, Z) or (G, I) plane from the point
  % \gamma_{\vartheta_{old}} to the stable manifold of q.
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
  %          * u_in(4:6) - Initial point of the manifold segment.
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
  pdim = data_in.pdim;

  %---------------%
  %     Input     %
  %---------------%
  % Initial point of the periodic orbit
  x0_seg = u_in(1 : xdim);
  % Final point of the periodic orbit
  x0_Wsq = u_in(xdim+1 : 2*xdim);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector of periodic orbit segment
  vec_PO = [x0_seg(1), x0_seg(3)];
  % Vector of stable manifold segment
  vec_W  = [x0_Wsq(1), x0_Wsq(3)];

  % Calculate difference
  vec_diff = vec_W - vec_PO;

  % Displacement amplitude
  A_perturb = norm(vec_diff);

  % Angle of displacement vector
  theta_perturb = atan2(vec_diff(2), vec_diff(1));
  % theta_perturb = mod(theta_perturb, 2*pi);
  theta_perturb = theta_perturb / (2 * pi);

  %----------------%
  %     Output     %
  %----------------%
  y_out = [A_perturb; theta_perturb];

end