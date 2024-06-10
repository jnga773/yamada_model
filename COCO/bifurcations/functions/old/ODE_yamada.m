function y_out = ODE_yamada(t_in, u_in, p_in)
  % ODE_YAMADA: Time-dependent encoding of the Yamada model set
  % of equations.

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Grab the state-space variables from u_in
  % Gain
  G = u_in(1, :);
  % Absorption
  Q = u_in(2, :);
  % Laser intensity
  I = u_in(3, :);

  % Grab the parameter-space variables from p_in
  % Decay time of gain
  gamma = p_in(1);
  % Pump current on the gain
  A = p_in(2);
  % (Relative) absorption
  B = p_in(3);
  a = p_in(4);
    
  % State array
  x = [G; Q; I];
  % Parameter array
  p = [gamma; A; B; a];

  %----------------%
  %     Output     %
  %----------------%
  % The system of equations
  y_out = yamada(x, p);

end