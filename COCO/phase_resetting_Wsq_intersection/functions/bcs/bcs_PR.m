function [data_in, y_out] = bcs_PR(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_PR(prob_in, data_in, u_in)
  %
  % Boundary conditions for the four segments of the phase-resetting problem:
  %                          x1(0) - x2(1) = 0 ,
  %                          x1(1) - x2(0) = 0 ,
  %                          e1 . F(x1(0)) = 0 .
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
  %            u_in(1:3)   - x(0) of segment 1,
  %            u_in(4:6)   - x(0) of segment 2,
  %            u_in(7:9)   - x(1) of segment 1,
  %            u_in(10:12) - x(1) of segment 2,
  %            u_in(13:17) - Parameters.
  %
  % Returns
  % -------
  % y_out : array of vectors
  %     An array containing the boundary conditions.
  % data_in : structure
  %     Function data structure to give dimensions of parameter and state
  %     space.

  % (defined in calc_PR_initial_conditions.m)
  % Original vector space dimensions
  xdim   = data_in.xdim;
  pdim   = data_in.pdim;
  % Parameter maps
  p_maps = data_in.p_maps;

  %============================================================================%
  %                              INPUT PARAMETERS                              %
  %============================================================================%
  %--------------------------------%
  %     Input: Initial Vectors     %
  %--------------------------------%
  % Segment 1 - x(0)
  x0_seg1       = u_in(1 : xdim);
  % Segment 2 - x(0)
  x0_seg2       = u_in(xdim+1 : 2*xdim);

  %------------------------------%
  %     Input: Final Vectors     %
  %------------------------------%
  % Segment 1 - x(1)
  x1_seg1       = u_in(2*xdim+1 : 3*xdim);
  % Segment 2 - x(1)
  x1_seg2       = u_in(3*xdim+1 : 4*xdim);
  %---------------------------%
  %     Input: Parameters     %
  %---------------------------%
  % Parameters
  parameters    = u_in(4*xdim+1 : end);

  % System parameters
  p_system     = parameters(1 : pdim);

  % Phase resetting parameters
  theta_gamma = parameters(p_maps.theta_gamma);

  %============================================================================%
  %                         BOUNDARY CONDITION ENCODING                        %
  %============================================================================%
  %---------------------------------%
  %     Segment 1 and Segment 2     %
  %---------------------------------%
  % Identity matrix
  ones_matrix = eye(xdim);
  % First component unit vector
  e1 = ones_matrix(1, :);

  % Boundary Conditions - Segments 1 and 2
  bcs_seg12_1   = x0_seg1 - x1_seg2;
  bcs_seg12_2   = x1_seg1 - x0_seg2;
  bcs_seg12_3   = e1 * yamada(x0_seg1, p_system);

  %============================================================================%
  %                                   OUTPUT                                   %
  %============================================================================%
  %----------------%
  %     Output     %
  %----------------%
  y_out = [bcs_seg12_1; cs_seg12_2; bcs_seg12_3;];

end
