function y_out = func_seg1(x_in, p_in)
  % y_out = func_seg1(u_in, p_in)
  %
  % COCO 'ode' toolbox encoding for the vector field corresponding to
  % segment 1 of the phase resetting curve.
  %
  % Segment 1 goes from gamma_0 to theta_new.
  %
  % Input
  % ----------
  % x_in : array, float
  %     State vector for the periodic orbit (x) and perpendicular
  %     vector (w).
  % p_in : array, float
  %     Array of parameter values
  %
  % Output
  % ----------
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
  x_vec        = x_in(1:xdim, :);

  % Perpendicular vectors
  w_vec        = x_in(xdim+1:2*xdim, :);

  % System parameters
  p_system     = p_in(1:pdim, :);

  % Phase resetting parameters
  % Period of the segment
  T             = p_in(pdim+1, :);
  % Integer for period
  % k             = p_in(pdim+2, :);
  k             = 1;
  % Stable Floquet eigenvalue
  % mu_s          = p_in(pdim+3, :);
  % Distance from pertured segment to \Gamma
  % eta           = p_in(pdim+4, :);
  % Phase where perturbation starts
  % theta_old     = p_in(pdim+5, :);
  % Phase where segment comes back to \Gamma
  theta_new     = p_in(pdim+6, :);
  % Angle of perturbation
  % theta_perturb = p_in(pdim+7, :);
  % Azimuthal angle of perturbation
  % phi_perturb   = p_in(pdim+8, :);
  % Size of perturbation
  % A_perturb     = p_in(pdim+9, :);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Calculate vector field
  vec_field = field(x_vec, p_system);
  
  % Save to array
  % vec_eqn = k .* theta_new .* vec_field;
  vec_eqn = k .* T .* theta_new .* vec_field;

  % Calculate adjoint equations
  % Jacobian at the zero-phase point
  J = field_DFDX(x_vec, p_system);

  % Cycle through each variable in x1 and calculate
  % adjoint equation components
  for i = 1 : numel(x_vec(1, :))
    % Transpose of Jacobian
    J_transpose(:, :, i) = transpose(J(:, :, i));

    % Calculate some things
    temp(:, :, i) = -theta_new(i) * J_transpose(:, :, i);

    % Save to array
    % adj_eqn(:, :, i) = temp(:, :, i) * w_vec(:, i);
    adj_eqn(:, :, i) = T(i) * temp(:, :, i) * w_vec(:, i);
  end

  %----------------%
  %     Output     %
  %----------------%
  % Vector field
  y_out(1:xdim, :) = vec_eqn(:, :);


  % Adjoint equation
  y_out(xdim+1:2*xdim, :) = adj_eqn(:, :);

end