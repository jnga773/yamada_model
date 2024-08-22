function plot_base_periodic_orbit(ax_in)
  % Plot the solution g
  
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Load data matrix
  load('./data_mat/initial_PO_manifold_q.mat', 'Wq_s', 'xbp', 'xpos', 'xneg', 'x0');

  %--------------%
  %     Plot     %
  %--------------%
  % Plotting colours
  colours = colororder();

  % Plot initial periodic orbit
  plot3(ax_in, xbp(:, 1), xbp(:, 2), xbp(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

  % Plot equilibrium points: x_{+}
  plot3(ax_in, xpos(1), xpos(2), xpos(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
        MarkerEdgeColor='b', DisplayName='$q$');

  % Plot equilibrium points: x_{-}
  plot3(ax_in, xneg(1), xneg(2), xneg(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
        MarkerEdgeColor='r', DisplayName='$p$');

  % Plot equilibrium points: x_{0}
  plot3(ax_in, x0(1), x0(2), x0(3), ...
        LineStyle='none', ...
        Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
        MarkerEdgeColor='r', DisplayName='$o$');
  
  % Plot stable manifold of q / x_{+}
  plot3(ax_in, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
        Color=colours(1, :), ...
        DisplayName='$W^{s}(p)$');

end