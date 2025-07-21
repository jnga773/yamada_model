function y_out = func_seg1(x_in, p_in)
  % y_out = func_seg1(u_in, p_in)
  %
  % COCO 'ode' toolbox encoding for the vector field corresponding to
  % segment 1 of the phase resetting curve.
  %
  % Segment 1 goes from gamma_0 to theta_new.
  %
  % Parameters
  % ----------
  % x_in : array, float
  %     State vector for the periodic orbit (x) and perpendicular
  %     vector (w).
  % p_in : array, float
  %     Array of parameter values
  %
  % Returns
  % -------
  % y_out : array, float
  %     Array of the vector field of the periodic orbit segment
  %     and the corresponding adjoint equation for the perpendicular
  %     vector.

  % Original vector field dimensions (CHANGE THESE)
  xdim = 3;
  pdim = 4;
  % Original vector field function
  field      = @yamada;
  field_DFDX = @yamada_DFDX;

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % State space variables
  x_vec        = x_in(1:xdim, :)

  % System parameters
  p_system     = p_in(1:pdim, :);

  % Phase resetting parameters
  theta_gamma = p_in(pdim+1, :);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Calculate vector field
  vec_field = field(x_vec, p_system);
  
  % Save to array
  vec_eqn = theta_gamma .* vec_field;

  %----------------%
  %     Output     %
  %----------------%
  % Vector field
  y_out(1:xdim, :) = vec_eqn(:, :);

end
