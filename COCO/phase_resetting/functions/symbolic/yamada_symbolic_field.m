function F_out = yamada_symbolic_field()
  % F_out = yamada_symbolic_field()
  %
  % Symbolic notation of the Yamada vector field in the transformed axes.
  
  % State-space and parameter variables
  syms G Q I
  syms gam A B a

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector field components
  % d/dt G
  F1 = gam * (A - G - (G * I));

  % d/dt Q
  F2 = gam * (B - Q - (a * Q * I));

  % d/dt QI
  F3 = I * (G - Q - 1);

  % Vector field
  F_vec = [F1; F2; F3];

  %----------------%
  %     Output     %
  %----------------%
  F_out = F_vec;

end