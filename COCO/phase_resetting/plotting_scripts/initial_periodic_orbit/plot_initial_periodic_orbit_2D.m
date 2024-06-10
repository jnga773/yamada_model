function plot_initial_periodic_orbit_2D(save_figure)
  % plot_initial_periodic_orbit_2D()
  %
  % Plots the initial periodic orbit from the 'coll' toolbox run.
  % along with the three stationary points ('q', 'p', and 'o'),
  % and the one-dimensional stable manifold of point 'q'. The data
  % is read from the file './data/initial_PO.mat'.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Load data matrix
  load('./data/initial_PO.mat');

  % Parameters
  A = p_PO(2);
  B = p_PO(3);

  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
  % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(3); fig.Name = 'Initial Periodic Orbits'; clf;
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
  plot(ax, xbp_PO(:, 1), xbp_PO(:, 2), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

  % Plot equilibrium points: x_{+}
  plot(ax, x_pos(1), x_pos(2), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
        MarkerEdgeColor='b', DisplayName='$q$');

  % Plot equilibrium points: x_{-}
  plot(ax, x_neg(1), x_neg(2), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
        MarkerEdgeColor='r', DisplayName='$p$');

  % Plot equilibrium points: x_{0}
  plot(ax, x_0(1), x_0(2), ...
        LineStyle='none', ...
        Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
        MarkerEdgeColor='r', DisplayName='$o$');

  % Legend
  legend(ax, Interpreter='latex');

  % Turn of axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % % X-Axis
  % ax.XAxis.TickDirection = 'in';
  % ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
  % ax.XAxis.MinorTick = 'on';
  % ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;
  % 
  % % Y-Axis
  % ax.YAxis.TickDirection = 'in';
  % ax.YAxis.TickValues = 0.0 : 2.0 : 8.0;
  % ax.YAxis.MinorTick = 'on';
  % ax.YAxis.MinorTickValues = 1.0 : 2.0 : 8.0;
  % 
  % % Z-Axis
  % ax.ZAxis.TickDirection = 'in';
  % ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
  % ax.ZAxis.MinorTick = 'on';
  % ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  % ax.XAxis.Limits = [0.0, 10.0];
  % ax.YAxis.Limits = [0.0, 8.0];
  % ax.ZAxis.Limits = [0.0, 18.0];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$G(t)$';
  ax.YAxis.Label.String = '$Q(t)$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Initial Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %----------------------%
  %      Save Figure     %
  %----------------------%
  if save_figure == true
    % Filename
    figname = 'initial_periodic_orbit_2D';
    % exportgraphics(fig, ['./images/', figname, '.png'], Resolution=800);
    exportgraphics(fig, ['./images/', figname, '.pdf'], ContentType='vector');
  end

end