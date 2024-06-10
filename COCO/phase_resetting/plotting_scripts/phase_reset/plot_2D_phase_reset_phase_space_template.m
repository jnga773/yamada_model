function plot_2D_phase_reset_phase_space_template(ax_in, run_in, label_in)
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
  
  % Read equilibrium point data
  sol_ep = ep_read_solution('singularity', run_in, label_in);

  % Equiliubrium point
  x_ep = sol_ep.x;

  % Read original periodic orbit solution
  sol_PO = coll_read_solution('initial_PO', 'run06_initial_periodic_orbit', 1);

  % Initial periodic orbit
  load('./data/initial_PO.mat');

  %-------------------%
  %     Plot Data     %
  %-------------------%
  % Default line colours
  colours = colororder();

  % Plot base solution
  plot(ax_in, xbp_PO(:, 1), xbp_PO(:, 3), ...
       LineStyle='-', Color=colours(3, :), ...
       DisplayName='$\Gamma$');

  % Plot segment 1
  plot(ax_in, x1(:, 1), x1(:, 3), Color=colours(1, :), ...
       DisplayName='Segment 1');

  % Plot segment 2
  plot(ax_in, x2(:, 1), x2(:, 3), Color=colours(2, :), ...
       DisplayName='Segment 2');

  % Plot segment 3
  plot(ax_in, x3(:, 1), x3(:, 3), Color=colours(3, :), ...
       DisplayName='Segment 3');

  % Plot segment 4
  plot(ax_in, x4(:, 1), x4(:, 3), Color=[colours(4, :), 0.25], ...
       DisplayName='Segment 4');

  % Plot start point of segment 4
  plot(ax_in, x4(1, 1), x4(1, 3), LineStyle='none', Marker='square', MarkerSize=10, ...
       MarkerFaceColor=colours(4, :), MarkerEdgecolor=colours(4, :), ...
       DisplayName='$x_{4}(0)$');
  % Plot end point of segment 4
  plot(ax_in, x4(end, 1), x4(end, 3), LineStyle='none', Marker='diamond', MarkerSize=10, ...
       MarkerFaceColor=colours(4, :), MarkerEdgecolor=colours(4, :), ...
       DisplayName='$x_{4}(1)$');

% Plot equilibrium points: x_{+}
  plot(ax_in, x_pos(1), x_pos(3), ...
       LineStyle='none', ...
       Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
       MarkerEdgeColor='b', DisplayName='$q$');

  % Plot equilibrium points: x_{-}
  plot(ax_in, x_neg(1), x_neg(3), ...
       LineStyle='none', ...
       Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
       MarkerEdgeColor='r', DisplayName='$p$');

  % Plot equilibrium points: x_{0}
  plot(ax_in, x_0(1), x_0(3), ...
       LineStyle='none', ...
       Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
       MarkerEdgeColor='r', DisplayName='$o$');
  
  % % Plot stable manifold of q / x_{+}
  % plot(ax_in, W_q_stable(:, 1), W_q_stable(:, 3), ...
  %      Color=colours(1, :), ...
  %      DisplayName='$W^{s}(p)$');


end