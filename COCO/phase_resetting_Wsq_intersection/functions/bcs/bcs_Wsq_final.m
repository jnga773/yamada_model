function [data_in, y_out] = bcs_Wsq_final(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_Wsq_final(prob_in, data_in, u_in)
  % 
  % COCO compatible encoding for the "initial" boundary conditions of the
  % trajectory segment of the stable manifold of x_pos. The segment starts
  %  near the equilibrium point, with boundary condition:
  %           x1(1) = x_pos + (epsilon * vs) .
  % Here, epsilon is the distances from of the trajectory to the equilbrium
  % points and vs is the stable eigenvector of the Jacobian at x_pos.
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
  %          * u_in(1:3) - The final point of the manifold trajectory segment (x1_Wsq),
  %          * u_in(4:6) - The x_pos equilibrium point (x_pos),
  %          * u_in(7:9) - Stable eigenvector (vs),
  %          * u_in(10)  - The epsilon spacings (eps1, eps2).
  %
  % Output
  % ----------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Not actually output here but you need to have it for COCO.

  % Original vector field dimensions
  xdim = data_in.xdim;

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Final vector of segment 1
  x1_Wsq  = u_in(1 : xdim);
  % x_pos equilibrium point
  xpos    = u_in(xdim+1 : 2*xdim);
  % Epsilon spacings and angle
  epsilon = u_in(end);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Eigenvector matrix of the Jacobian of x_pos
  eigvec = data_in.ep_X;

  % Get stable eigenvector
  vs = eigvec(:, 3);

  %----------------%
  %     Output     %
  %----------------%
  % Boundary conditions
  y_out = x1_Wsq - (xpos + (epsilon * vs));

end