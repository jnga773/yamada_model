function plot_initial_periodic_orbit_2D_quasi(save_figure)
  % plot_initial_periodic_orbit_2D_quasi()
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

  % Transform data
  xbp_PO_2D = [(A * xbp_PO(:, 1) / ((A ^ 2) + (B ^ 2))) + (B * xbp_PO(:, 2) / ((A ^ 2) + (B ^ 2))), xbp_PO(:, 3)];
  x_neg_2D  = [(A * x_neg(1) / ((A ^ 2) + (B ^ 2))) + (B * x_neg(2) / ((A ^ 2) + (B ^ 2))), x_neg(3)];
  x_pos_2D  = [(A * x_pos(1) / ((A ^ 2) + (B ^ 2))) + (B * x_pos(2) / ((A ^ 2) + (B ^ 2))), x_pos(3)];
  x_0_2D    = [(A * x_0(1) / ((A ^ 2) + (B ^ 2))) + (B * x_0(2) / ((A ^ 2) + (B ^ 2))), x_0(3)];

  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
  % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(2); fig.Name = 'Initial Periodic Orbits'; clf;
  fig.Units = 'inches'; fig.Position = [3, 3, 8, 8]; fig.PaperSize = [8, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plot initial periodic orbit
  plot(ax, xbp_PO_2D(:, 1), xbp_PO_2D(:, 2), ...
       LineStyle='-', Color=colours(3, :), ...
       DisplayName='$\Gamma$');

  % Plot equilibrium points: x_{+}
  plot(ax, x_pos_2D(1), x_pos_2D(2), ...
       LineStyle='none', ...
       Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
       MarkerEdgeColor='b', DisplayName='$q$');

  % Plot equilibrium points: x_{-}
  plot(ax, x_neg_2D(1), x_neg_2D(2), ...
       LineStyle='none', ...
       Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
       MarkerEdgeColor='r', DisplayName='$p$');

  % Plot equilibrium points: x_{0}
  plot(ax, x_0_2D(1), x_0_2D(2), ...
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
  % X-Axis
  % ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
  % ax.XAxis.MinorTick = 'on';
  % ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;
  
  % Y-Axis
  % ax.YAxis.TickValues = 0.0 : 2.0 : 8.0;
  % ax.YAxis.MinorTick = 'on';
  % ax.YAxis.MinorTickValues = 1.0 : 2.0 : 8.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  % ax.XAxis.Limits = [0.0, 10.0];
  % ax.YAxis.Limits = [0.0, 8.0];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$\tilde{G}(t)$';
  ax.YAxis.Label.String = '$I(t)$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Initial Periodic Orbit (2D Transformed Axes)';

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
    figname = 'initial_periodic_orbit_2D_quasi';
    % exportgraphics(fig, ['./images/', figname, '.png'], Resolution=800);
    exportgraphics(fig, ['./images/', figname, '.pdf'], ContentType='vector');
  end

end

function x_out = transform_axes(x_in, A_in, B_in)
  % transforms the G coordinate

  % G component
  G1 = (A_in * x_in(:, 1)) / (A_in ^ 2) + (B_in ^ 2);
  % Q compnent
  G2 = (B_in * x_in(:, 2)) / (A_in ^ 2) + (B_in ^ 2);

  % Transformed G
  G = G1 + G2;

  % Output
  x_out = [G, x_in(:, 3)];

end