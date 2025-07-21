function plot_base_periodic_orbit(ax_in)
  % Plot the solution g
  
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Load data matrix
  % load('./data_mat/PO_and_manifolds.mat', 'xbp_PO', ...
  %      'Ws_q', 'Ws_PO1', 'Ws_PO2', ...
  %      'xpos', 'xneg', 'x0');
  load('./data_mat/initial_PO.mat');

  % Plotting colours
  colours = colororder();

  %------------------------------%
  %     Plot: Periodic Orbit     %
  %------------------------------%
  % Plot initial periodic orbit
  plot3(ax_in, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

  %---------------------------------%
  %     Plot: Stationary Points     %
  %---------------------------------%
  % Plot equilibrium points: x_{+}
  plot3(ax_in, xpos(1), xpos(2), xpos(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
        MarkerEdgeColor='b', DisplayName='$q$');

  % % Plot equilibrium points: x_{-}
  % plot3(ax_in, xneg(1), xneg(2), xneg(3), ...
  %       LineStyle='none', ...
  %       Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
  %       MarkerEdgeColor='r', DisplayName='$p$');
  % 
  % % Plot equilibrium points: x_{0}
  % plot3(ax_in, x0(1), x0(2), x0(3), ...
  %       LineStyle='none', ...
  %       Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
  %       MarkerEdgeColor='r', DisplayName='$o$');

  %--------------------------------%
  %     Plot: Stable Manifolds     %
  %--------------------------------%
  % Plot stable manifold of q / x_{+}
  plot3(ax_in, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), Color=colours(1, :), ...
        DisplayName='$W^{s}(q)$');
  
end