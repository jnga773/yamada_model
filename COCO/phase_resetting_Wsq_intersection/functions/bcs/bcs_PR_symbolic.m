function bcs_coco_out = bcs_PR_symbolic()
  % bcs_coco_out = bcs_PR_symbolic()
  %
  % Boundary conditions for the four segments of the phase-resetting problem:
  %                          x1(0) - x2(1) = 0 ,
  %                          x1(1) - x2(0) = 0 ,
  %                          e1 . F(x1(0)) = 0 .
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %            u_in(1:3)   - x(0) of segment 1,
  %            u_in(4:6)   - x(0) of segment 2,
  %            u_in(7:9)   - x(1) of segment 1,
  %            u_in(10:12) - x(1) of segment 2,
  %            u_in(13:17) - Parameters.
  %
  % Returns
  % -------
  % bcs_coco_out : cell of function handles
  %     List of CoCo-ified symbolic functions for the boundary conditions
  %     Jacobian, and Hessian.

  % State-space dimension
  xdim = 3;

  %============================================================================%
  %                              INPUT PARAMETERS                              %
  %============================================================================%
  %--------------------------------%
  %     Input: Initial Vectors     %
  %--------------------------------%
  % Segment 1 - x(0)
  x0_seg1 = sym('x0_seg1', [xdim, 1]);
  % Segment 2 - x(0)
  x0_seg2 = sym('x0_seg2', [xdim, 1]);

  %------------------------------%
  %     Input: Final Vectors     %
  %------------------------------%
  % Segment 1 - x(1)
  x1_seg1 = sym('x1_seg1', [xdim, 1]);
  % Segment 2 - x(1)
  x1_seg2 = sym('x1_seg2', [xdim, 1]);

  %---------------------------%
  %     Input: Parameters     %
  %---------------------------%
  % System parameters
  syms gam A B a
  p_sys = [gam; A; B; a];

  % Phase resetting parameters
  syms theta_gamma
  p_PR = theta_gamma;

  %============================================================================%
  %                         BOUNDARY CONDITION ENCODING                        %
  %============================================================================%
  %---------------------------------%
  %     Segment 1 and Segment 2     %
  %---------------------------------%
  % Vector field
  F_vec = yamada_symbolic_field(x0_seg1, p_sys);

  % Boundary Conditions - Segments 1 and 2
  bcs_seg12_1   = x0_seg1 - x1_seg2;
  bcs_seg12_2   = x0_seg2 - x1_seg1;
  bcs_seg12_3   = F_vec(1);

  %============================================================================%
  %                                   OUTPUT                                   %
  %============================================================================%
  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Combined vector
  uvec = [x0_seg1; x0_seg2; x1_seg1; x1_seg2;
          p_sys; p_PR];

  % Boundary conditions vector
  bcs =  [bcs_seg12_1;  bcs_seg12_2; bcs_seg12_3];

  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_PR';

  % COCO Function encoding
  bcs_coco = sco_sym2funcs(bcs, {uvec}, {'u'}, 'filename', filename_out);

  % Function to "CoCo-ify" function outputs: [data_in, y_out] = f(prob_in, data_in, u_in)
  cocoify = @(func_in) @(prob_in, data_in, u_in) deal(data_in, func_in(u_in));

  % List of functions
  func_list = {cocoify(bcs_coco('')), cocoify(bcs_coco('u')), cocoify(bcs_coco({'u', 'u'}))};

  %----------------%
  %     Output     %
  %----------------%
  bcs_coco_out = func_list;

end
