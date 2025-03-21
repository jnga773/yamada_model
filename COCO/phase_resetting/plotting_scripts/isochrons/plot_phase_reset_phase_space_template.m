function plot_phase_reset_phase_space_template(ax_in, run_in, label_in)
  % plot_phase_reset_phase_space_template(ax_in, run_in, label_in)
  %
  % Plots the state space solution from solution "label_in" of "run_in"
  % to axis "ax_in".

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read segment data
  sol1 = coll_read_solution('seg1', run_in, label_in);
  sol2 = coll_read_solution('seg2', run_in, label_in);
  sol3 = coll_read_solution('seg3', run_in, label_in);
  sol4 = coll_read_solution('seg4', run_in, label_in);

  % State space data
  x1   = sol1.xbp;
  x2   = sol2.xbp;
  x3   = sol3.xbp;
  x4   = sol4.xbp;

  %-------------------%
  %     Plot Data     %
  %-------------------%
  % Default line colours
  colours = colororder();

  % Plot base solution
  plot_base_periodic_orbit(ax_in, false);

  % Plot segment 1
  plot3(ax_in, x1(:, 1), x1(:, 2), x1(:, 3), Color=colours(1, :), ...
        DisplayName='Segment 1');

  % Plot segment 2
  plot3(ax_in, x2(:, 1), x2(:, 2), x2(:, 3), Color=colours(2, :), ...
        DisplayName='Segment 2');

  % Plot segment 3
  plot3(ax_in, x3(:, 1), x3(:, 2), x3(:, 3), Color=colours(3, :), ...
        DisplayName='Segment 3');

  % Plot segment 4
  plot3(ax_in, x4(:, 1), x4(:, 2), x4(:, 3), Color=[colours(4, :), 0.7], ...
        DisplayName='Segment 4');

  % Plot start point of segment 4
  plot3(ax_in, x4(1, 1), x4(1, 2), x4(1, 3), LineStyle='none', Marker='square', MarkerSize=10, ...
        MarkerFaceColor=colours(4, :), MarkerEdgecolor=colours(4, :), ...
        DisplayName='$x_{4}(0)$');
  % Plot end point of segment 4
  plot3(ax_in, x4(end, 1), x4(end, 2), x4(end, 3), LineStyle='none', Marker='diamond', MarkerSize=10, ...
        MarkerFaceColor=colours(4, :), MarkerEdgecolor=colours(4, :), ...
        DisplayName='$x_{4}(1)$');

end
