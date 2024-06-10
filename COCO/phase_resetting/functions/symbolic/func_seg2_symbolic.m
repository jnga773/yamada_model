function F_coco_out = func_seg2_symbolic()
  % F_coco_out = func_seg2_symbolic()
  %
  % Creates a CoCo-compatible function encoding for the first
  % segment of the phase-resetting problem.
  %
  % Segment 2 goes from theta_new to gamma_0.

  % State space dimension
  xdim = 3;

  %---------------%
  %     Input     %
  %---------------%
  % State-space variables
  syms G Q I
  xvec = [G; Q; I];

  % Adjoint equation variables
  wvec = sym('w', [xdim, 1]);

  % System parameters
  syms gam A B a
  p_sys = [gam; A; B; a];

  % Phase resetting parameters
  syms T k mu_s eta
  syms theta_old theta_new
  syms theta_perturb phi_perturb A_perturb
  p_PR = [T; k; mu_s; eta;
          theta_old; theta_new;
          theta_perturb; phi_perturb; A_perturb];

  % Total vectors
  uvec = [xvec; wvec];
  pvec = [p_sys; p_PR];

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector field
  F_vec = yamada_symbolic_field();

  % Vector field equation
  % vec_eqn = (1 - theta_new) * F_vec;
  vec_eqn = T * (1 - theta_new) * F_vec;

  % Calculate tranpose of Jacobian at point xvec
  J_T = transpose(jacobian(F_vec, xvec));

  % Adjoint equation
  % adj_eqn = -(1 - theta_new) * J_T * wvec;
  adj_eqn = -T * (1 - theta_new) * J_T * wvec;

  % Total equation
  F_seg = [vec_eqn; adj_eqn];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symbolic/F_seg2';

  % COCO Function encoding
  F_coco = sco_sym2funcs(F_seg, {uvec, pvec}, {'x', 'p'}, 'filename', filename_out);

  % List of functions
  func_list = {F_coco(''), ...
               F_coco('x'), F_coco('p'), ...
               F_coco({'x', 'x'}), F_coco({'x', 'p'}), F_coco({'p', 'p'})};

  %----------------%
  %     Output     %
  %----------------%
  F_coco_out = func_list;

end