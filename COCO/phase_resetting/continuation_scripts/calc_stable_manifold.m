function x_out = calc_stable_manifold(run_in, label_in)
  % x_out = calc_stable_manifold(run_in, label_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read EP solution
  sol_EP   = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium point
  x_pos = sol_EP.x;
  % Parameters
  p     = sol_EP.p;
  
  %------------------------------%
  %     Calculate EigenStuff     %
  %------------------------------%
  % Jacobian
  J_stable = yamada_DFDX(x_pos, p);

  % Calculate eigenvalues and eigenvectors
  [eigvec, eigval] = eig(J_stable);

  % Indices for stable eigenvectors (eigval < 0)
  % stable_index = find(diag(eigval) < 0);
  stable_index = 3;

  % Stable eigenvector
  vec_s = eigvec(:, stable_index);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps1 = -0.01;
  % Time span
  % t_span1 = [0.0, -17.0];
  t_span1 = -16.5 : 0.01 : 0.0;
  t_span1 = flip(t_span1);

  % Initial vector
  x_init1 = x_pos + (eps1 * vec_s);

  % Integrate using ode45
  [~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span1, x_init1);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps2 = 0.01;
  % Time span
  % t_span2 = [0.0, -25.0];
  t_span2 = -25.0 : 0.01 : 0.0;
  t_span2 = flip(t_span2);

  % Initial vector
  x_init2 = x_pos + (eps2 * vec_s);

  % Integrate using ode45
  [~, W2] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span2, x_init2);

  %----------------%
  %     Output     %
  %----------------%
  x_out = [flip(W2); W1];

end