function plot_initial_periodic_orbit_COLL()
  % plot_initial_periodic_orbit_COLL()
  %
  % Plots the initial periodic orbit from the 'coll' toolbox run.
  % along with the three stationary points ('q', 'p', and 'o'),
  % and the one-dimensional stable manifold of point 'q'. The data
  % is read from the file './data/initial_PO.mat'.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Load data matrix
  load('./data_mat/initial_PO.mat');

  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
  % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(1); fig.Name = 'Initial Periodic Orbits'; clf;
  fig.Units = 'inches'; fig.Position = [3, 3, 8, 8]; fig.PaperSize = [8, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plotting colours
  colours = colororder();

  % Plot initial periodic orbit
  plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

  % Plot equilibrium points: x_{+}
  plot3(ax, x_pos(1), x_pos(2), x_pos(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
        MarkerEdgeColor='b', DisplayName='$q$');

  % Plot equilibrium points: x_{-}
  plot3(ax, x_neg(1), x_neg(2), x_neg(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
        MarkerEdgeColor='r', DisplayName='$p$');

  % Plot equilibrium points: x_{0}
  plot3(ax, x_0(1), x_0(2), x_0(3), ...
        LineStyle='none', ...
        Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
        MarkerEdgeColor='r', DisplayName='$o$');
  
  % Plot stable manifold of q / x_{+}
  plot3(ax, W_q_stable(:, 1), W_q_stable(:, 2), W_q_stable(:, 3), ...
        Color=colours(1, :), ...
        DisplayName='$W^{s}(p)$');

  % Legend
  legend(ax, Interpreter='latex');

  % Turn of axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickDirection = 'in';
  ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;
  
  % Y-Axis
  ax.YAxis.TickDirection = 'in';
  ax.YAxis.TickValues = 0.0 : 2.0 : 8.0;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = 1.0 : 2.0 : 8.0;
  
  % Z-Axis
  ax.ZAxis.TickDirection = 'in';
  % ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
  % ax.ZAxis.MinorTick = 'on';
  % ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 10.0];
  ax.YAxis.Limits = [0.0, 8.0];
  ax.ZAxis.Limits = [0.0, ceil(max(xbp_PO(:, 3)))];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$G(t)$';
  ax.YAxis.Label.String = '$Q(t)$';
  ax.ZAxis.Label.String = '$I(t)$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Initial Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % 3D plot view
  view(45, 15.0);

  %----------------------%
  %      Save Figure     %
  %----------------------%
  % Filename
  figname = 'initial_periodic_orbit_3D';
  exportgraphics(fig, ['./images/', figname, '.pdf'], ContentType='vector');

end