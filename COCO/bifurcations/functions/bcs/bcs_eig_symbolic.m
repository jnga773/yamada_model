function bcs_coco_out = bcs_eig_symbolic()
  % bcs_coco_out = bcs_eig_symbolic()
  %
  % COCO compatible encoding for the boundary conditions of the eigenvalues and
  % eigenvectors of the Jacobian evaluated at the stationary point x_neg
  %
  % For the hardcoded version, and the actual functions that
  % will be coco_add_func call will include the following
  % u-vector components:
  %          u_in(1:3)  - Equilibrium point
  %          u_in(4:7)  - System parameters
  %          u_in(8:10) - The eigenvector
  %          u_in(11)   - The eigenvalue
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
  % Initial point of the periodic orbit
  x_neg = sym('x_neg', [xdim, 1]);

  % System parameters
  syms gam A B a
  p_sys = [gam; A; B; a];

  % Eigenvector
  eig_vec = sym('eigvec', [xdim, 1]);
  % Eigenvalue
  syms eig_val

  % Combined vector
  uvec = [x_neg;
          p_sys;
          eig_vec;
          eig_val];

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector field
  F_vec = yamada_symbolic_field(x_neg, p_sys);
  % Jacobian
  J     = jacobian(F_vec, x_neg);

  %---------------------------------------%
  %     Calculate Boundary Conditions     %
  %---------------------------------------%
  % Eigenvalue equations
  bcs1 = (J * eig_vec) - (eig_val * eig_vec)

  % Unit vector equations
  bcs2 = (eig_vec' * eig_vec) - 1

  % Boundary condition vector
  bcs = [bcs1; bcs2];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_bcs_eig';

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