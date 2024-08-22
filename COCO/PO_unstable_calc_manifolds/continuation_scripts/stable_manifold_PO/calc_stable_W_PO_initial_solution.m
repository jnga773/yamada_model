function data_out = calc_stable_W_PO_initial_solution(filename_in)
  % data_out = calc_stable_W_PO_initial_solution(filename_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read initial solution stuff from filename_in
  load(filename_in);

  %--------------------%
  %     Parameters     %
  %--------------------%
  % Period
  Tu    = T_u;
  % Angle
  theta = 1.0;

  % Extended parameter names
  pnames_u = {pnames{:}, 'Tu', 'theta'};

  %-------------------------------------%
  %     Setup Stable Manifold Stuff     %
  %-------------------------------------%
  % Initial distances from the equilibria, along the tangent spaces of the
  % unstable and stable manifolds, to the initial points along the corresponding
  % trajectory segments.
  eps = 0.01;

  % Initial time
  t0 = 0;

  % Initial state vector
  % x_init_1 = xbp_end_u' + (eps * vec_s);
  % x_init_2 = xbp_end_u' - (eps + vec_s);
  x_init_1 = xbp_u(1, :)';
  x_init_2 = xbp_u(1, :)';

  %----------------%
  %     Output     %
  %----------------%
  % Read periodic orbit solution
  data_out = load(filename_in);

  data_out.x_init_1 = x_init_1';
  data_out.x_init_2 = x_init_2';
  data_out.t0       = t0;

  data_out.lam_s    = lam_s;
  data_out.vec_s    = vec_s;

  data_out.eps      = eps;
  data_out.pnames_u = pnames_u;
  data_out.pu       = [p; Tu; theta];
  data_out.theta    = theta;
  data_out.tbp_u    = tbp_u / Tu;

end