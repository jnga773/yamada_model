function [data_in, y_out] = bcs_final(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_final(prob_in, data_in, u_in)
  %
  % Boundary conditions of the two trajectory segments. Both segments end on
  % the plane \Sigma, which we define by a point and a normal vector in
  % lins_method_setup().
  % The unstable manifold is solved in forwards time, with
  %           x_u(T1) \in \Sigma.
  % The stable manifold is solved in reverse time, with
  %           x_s(0) \in \Sigma.
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
  %          * u_in(1:3) - The final point of the unstable manifold (x1_unstable),
  %          * u_in(4:6) - The initial point of the stable manifold (x0_stable),
  %          * u_in(7:9) - The x_pos non-trivial equilibrium point in the
  %                        centre of the homoclinic; used to fix the z-position
  %                        of the \Sigma plane (x_pos).
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
  % Final vector of the unstable manifold
  x1_unstable = u_in(1 : xdim);
  % Initial vector of the stable manifold
  x0_stable   = u_in(xdim+1 : 2*xdim);
  % x_pos equilibrium point for \Sigma plane
  x_pos      = u_in(2*xdim+1 : end);

  %------------------------------------%
  %     Read Points from "data_in"     %
  %------------------------------------%
  % Point on plane \Sigma
  pt0 = data_in.pt0;
  pt0(xdim) = x_pos(xdim);

  % Normal vector of plane \Sigma
  normal = data_in.normal;

  %----------------%
  %     Output     %
  %----------------%
  % Boundary conditions
  y_out = [normal * (x1_unstable - pt0);
           normal * (x0_stable   - pt0)];

end