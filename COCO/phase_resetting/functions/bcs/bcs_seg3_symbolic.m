function bcs_coco_out = bcs_seg3_symbolic()
  % bcs_coco_out = bcs_seg3_symbolic()
  %
  % Boundary conditions for segment three of the phase reset
  % segments:
  %                        x3(1) - x1(0) = 0 .
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %           u_in(1:3)   - x(0) of segment 1,
  %           u_in(4:6)   - x(1) of segment 3.
  %
  % Output
  % ----------
  % bcs_coco_out : cell of function handles
  %     List of CoCo-ified symbolic functions for the boundary conditions
  %     Jacobian, and Hessian.

  % State-space dimension
  xdim = 3;

  %---------------%
  %     Input     %
  %---------------%
  % Segment 1 - x(0)
  x0_seg1 = sym('x0_seg1', [xdim, 1]);
  
  % Segment 3 - x(1)
  x1_seg3 = sym('x1_seg3', [xdim, 1]);

  % Combined vector
  uvec = [x0_seg1; x1_seg3];

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Boundary Conditions - Segment 3
  bcs_seg3 = x1_seg3 - x0_seg1;

  % Boundary condition vector
  bcs = bcs_seg3;

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_seg3';

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
