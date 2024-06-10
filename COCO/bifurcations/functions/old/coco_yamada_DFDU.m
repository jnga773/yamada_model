function [data_out, J_out] = COCO_yamada_DFDU(prob_in, data_out, u_in)
  % COCO_YAMADA_DFDU: COCO-compatible encoding of the Jacobian of
  % the Yamada model set of equations.

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

  %----------------%
  %     Output     %
  %----------------%
  % The Jacobian matrix of the state-space
  J_out = zeros(3, 7, numel(G));
  
  % State partial derivatives
  J_out(1, 1, :) = -gamma .* (1 + I);
  J_out(1, 2, :) = 0;
  J_out(1, 3, :) = -gamma .* G;

  J_out(2, 1, :) = 0;
  J_out(2, 2, :) = -gamma .* (1 + (a .* I));
  J_out(2, 3, :) = -gamma .* a .* Q;

  J_out(3, 1, :) = I;
  J_out(3, 2, :) = -I;
  J_out(3, 3, :) = G - Q - 1;

  % Parameter partial derivatives
  J_out(1, 4, :) = A - G - (G .* I);
  J_out(1, 5, :) = gamma;
  J_out(1, 6, :) = 0;
  J_out(1, 7, :) = 0;

  J_out(2, 4, :) = B - Q - (a .* Q .* I);
  J_out(2, 5, :) = 0;
  J_out(2, 6, :) = gamma;
  J_out(2, 7, :) = -gamma .* Q .* I;

  J_out(3, 4, :) = 0;
  J_out(3, 5, :) = 0;
  J_out(3, 6, :) = 0;
  J_out(3, 7, :) = 0;

end