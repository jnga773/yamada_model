function [data_out, y_out] = COCO_yamada(prob_in, data_out, u_in)
  % COCO_YAMADA: COCO-compatible encoding of the Yamada model
  % set of equations.

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
  gamma = u_in(4, :);
  % Pump current on the gain
  A = u_in(5, :);
  % (Relative) absorption
  B = u_in(6, :);
  a = u_in(7, :);
    
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