function F_out = yamada_symbolic_field()
  % F_out = yamada_symbolic_field()
  %
  % Symbolic notation of the Yamada vector field.
  
  % State-space and parameter variables
  % syms x1 x2 x3
  syms G Q I
  syms gam A B a

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  F1 = gam * (A - G - (G * I));
  F2 = gam * (B - Q - (a * Q * I));
  F3 = I * (G - Q - 1.0);

  % % Substitute (x1, x2, x3) for (G, Q, I)
  % F1 = subs(F1, [G, Q, I], [x1, x2, x3]);
  % F2 = subs(F2, [G, Q, I], [x1, x2, x3]);
  % F3 = subs(F3, [G, Q, I], [x1, x2, x3]);

  % Vector field
  F_vec = [F1; F2; F3];

  %----------------%
  %     Output     %
  %----------------%
  F_out = F_vec;

end