function F_out = yamada(x_in, p_in)
  % Symbolic encoding of the Yamada vector field
  
  % Get state vector components
  G = x_in(1);
  Q = x_in(2);
  I = x_in(3);

  % Get parameter components
  gamma = p_in(1);
  A     = p_in(2);
  B     = p_in(3);
  a     = p_in(4);

  % Vector field components
  F1 = gamma * (A - G - (G * I));
  F2 = gamma * (B - Q - (a * Q * I));
  F3 = I * (G - Q - 1);

  % Output
  F_out = [F1; F2; F3];

end
