function theta_new_out = fix_PTC_data(theta_new_in)
  % theta_new_out = fix_PTC_data(theta_new_in)
  %
  % If theta_new > 1, subtracts 1 until 0 <= theta_new <= 1.

  % Empty array
  theta_new_out = zeros(length(theta_new_in));

  % Cycle through data
  for i = 1 : length(theta_new_in)
    % Read data point
    theta_new = theta_new_in(i);

    % Minus the floor(theta_new)
    theta_new = theta_new - floor(theta_new);
    
    % Update array
    theta_new_out(i) = theta_new;
  end

end