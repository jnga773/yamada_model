function plot_homoclinic_hyperplane_base(ax_in, run_in, label_in, data_in)
  % PLOT_HOMOCLINIC_HYPERPLANE_BASE: Plots the unstable and stable
  % eigenvectors and the Sigma plane. This will be called in
  % plot_homoclinic_manifolds().

  %--------------------------------------------------%
  %     Read Data: Current Parameters and Points     %
  %--------------------------------------------------%
  [sol_neg, ~] = ep_read_solution('xneg', run_in, label_in);
  [sol_pos, ~] = ep_read_solution('xpos', run_in, label_in);

  % Parameters
  parameters = sol_neg.p;

  % Equilibrium points
  x0_neg = sol_neg.x;
  x0_pos = sol_pos.x;

  % Calculate non-trivial steady states and stable and unstable eigenvectors
  [vu, vs1, vs2] = unstable_stable_eigenvectors(x0_neg, parameters);

  % Scale for vectors
  scale = 0.1;

  plot_vu  = [
              [x0_neg(1), x0_neg(1) + (scale * vu(1))];
              [x0_neg(2), x0_neg(2) + (scale * vu(2))];
              [x0_neg(3), x0_neg(3) + (scale * vu(3))]
              ];
  plot_vs1 = [
              [x0_neg(1), x0_neg(1) + (scale * vs1(1))];
              [x0_neg(2), x0_neg(2) + (scale * vs1(2))];
              [x0_neg(3), x0_neg(3) + (scale * vs1(3))]
              ];
  plot_vs2 = [
              [x0_neg(1), x0_neg(1) + (scale * vs2(1))];
              [x0_neg(2), x0_neg(2) + (scale * vs2(2))];
              [x0_neg(3), x0_neg(3) + (scale * vs2(3))]
              ];

  %-------------------%
  %     Plot Data     %
  %-------------------%
  % Plot approximate homoclinic (high-period periodic orbit)
  plot3(ax_in, data_in.xbp_approx(:, 1), data_in.xbp_approx(:, 2), data_in.xbp_approx(:, 3), ...
        LineStyle='--', Color='Black', ...
        DisplayName='Approximate Homoclinic');

  % Plot manifold
  x = [-10, 10, 10, -10];
  y = [-10, -10, 10, 10];
  z = [x0_pos(3), x0_pos(3), x0_pos(3), x0_pos(3)];
  fill3(ax_in, x, y, z, 'blue', FaceColor='blue', FaceAlpha=0.25, ...
        DisplayName='$\Sigma$');

  % Plot dot for stationary point
  plot3(ax_in, x0_neg(1), x0_neg(2), x0_neg(3), Marker='o', LineStyle='none', ...
        MarkerSize=8, MarkerFaceColor='red', DisplayName='Equilibrium Point');

  % Plot vectors
  plot3(ax_in, plot_vu(1, :), plot_vu(2, :), plot_vu(3, :), '->', ...
        Color='blue', DisplayName='$\vec{v}^{u}$');
  plot3(ax_in, plot_vs1(1, :), plot_vs1(2, :), plot_vs1(3, :), '->', ...
        Color='red', DisplayName='$\vec{v}^{s}_{1}$');
  plot3(ax_in, plot_vs2(1, :), plot_vs2(2, :), plot_vs2(3, :), '->', ...
        Color='green', DisplayName='$\vec{v}^{s}_{2}$');

end