function [data_in, y_out] = bcs_isochron(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_isochron(prob_in, data_in, u_in)
  %
  % Boundary conditions for the isochron, that is:
  %            \theta_old - \theta_new = 0 .
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
  %          * u_in(1:3) - x(0) of segment 4
  %
  % Output
  % ----------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Function data structure to give dimensions of parameter and state
  %     space.

  % (defined in calc_PR_initial_conditions.m)
  % Original vector space dimensions
  xdim   = data_in.xdim;
  % pdim   = data_in.pdim;
  % % Parameter maps
  % p_maps = data_in.p_maps;

  %---------------%
  %     Input     %
  %---------------%
  % Segment 4 - x(0)
  x0_seg4 = u_in(1 : xdim);

  %----------------%
  %     Output     %
  %----------------%
  y_out = x0_seg4;

end