function [data_in, y_out] = bcs_initial(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_initial(prob_in, data_in, u_in)
  % 
  % COCO compatible encoding for the "initial" boundary conditions of the two
  % trajectory segments. Both segments start near the equilibrium point, x0.
  % The unstable manifold is solved in forwards time, with
  %           x_u(0) = x0 + eps1 * vu.
  %
  % The stable manifold is solved in reverse time, with
  %           x_s(T2) = x0 + eps2(vs1 * cos(theta) + vs2 * sin(theta)) .
  %
  % Here, eps1 and eps2 are the distances from of the trajectories to the equilbrium
  % points, vu is the unstable eigenvector of the Jacobian, vs1 and vs2 are the
  % stable eigenvectors of the Jacobian, and theta is the angle between the two.
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
  %          * u_in(1:3)   - The initial point of the unstable manifold (x0_unstable),
  %          * u_in(4:6)   - The final point of the stable manifold (x1_stable),
  %          * u_in(7:9)   - The x_neg equilibrium point (x_neg),
  %          * u_in(10:12) - Unstable eigenvector,
  %          * u_in(13:15) - Stable eigenvector 1,
  %          * u_in(16:18) - Stable eigenvector 2,
  %          * u_in(19:21) - The epsilon spacings and angle (eps).
  %
  % Output
  % ----------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Not actually output here but you need to have it for COCO.

  % State space dimension
  xdim = data_in.xdim;
  % Parameter space dimension
  % pdim = data_in.pdim;

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Initial vector of the unstable manifold
  x0_unstable = u_in(1 : xdim);
  % Final vector of the stable manifold
  x1_stable   = u_in(xdim+1 : 2*xdim);

  % x_neg equilibrium point
  x_neg      = u_in(2*xdim+1 : 3*xdim);
  % Unstable eigenvector
  vu          = u_in(3*xdim+1 : 4*xdim);
  % Stable eigenvectors
  vs1         = u_in(4*xdim+1 : 5*xdim);
  vs2         = u_in(5*xdim+1 : 6*xdim);

  % Epsilon spacings and angle
  eps1        = u_in(end-2);
  eps2        = u_in(end-1);
  theta       = u_in(end);

  %---------------------------------------%
  %     Calculate Boundary Conditions     %
  %---------------------------------------%
  % Unstable boundary condition
  x_init_u = x_neg + (eps1 * vu);

  % Stable boundary condition
  x_final_s = x_neg + eps2 * ((vs1 * cos(theta)) + (vs2 * sin(theta)));

  %----------------%
  %     Output     %
  %----------------%
  % Boundary conditions
  y_out = [x0_unstable - x_init_u ;
           x1_stable   - x_final_s];

end