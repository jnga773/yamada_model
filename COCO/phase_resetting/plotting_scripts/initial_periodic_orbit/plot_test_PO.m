function plot_test_PO(run_in, label_in)
  % Plot the solution g

  % Read da solution g
  sol = coll_read_solution('initial_PO', run_in, label_in);
  x_plot = sol.xbp;

  % Read da equilibrium points g
  sol1 = ep_read_solution('xpos', run_in, label_in);
  sol2 = ep_read_solution('xneg', run_in, label_in);
  sol3 = ep_read_solution('x0', run_in, label_in);

  xpos = sol1.x;
  xneg = sol2.x;
  x0   = sol3.x;

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
  
  % Plot base solution
  % plot_base_periodic_orbit(ax);

  % Plot continuation solution
  plot3(ax, x_plot(:, 1), x_plot(:, 2), x_plot(:, 3), LineStyle='-', Color=colours(1, :), ...
        DisplayName='Solution');

  % Plot equilibrium points
  plot3(ax, xpos(1), xpos(2), xpos(3), LineStyle='none', Marker='*', MarkerFaceColor='b', ...
        MarkerEdgeColor='b', DisplayName='$x_{+}^{*}$');
  plot3(ax, xneg(1), xneg(2), xneg(3), LineStyle='none', Marker='*', MarkerFaceColor='r', ...
        MarkerEdgeColor='r', DisplayName='$x_{-}^{*}$');
  plot3(ax, x0(1), x0(2), x0(3), LineStyle='none', Marker='o', MarkerFaceColor='r', ...
        MarkerEdgeColor='r', DisplayName='$x_{0}^{*}$');

  legend(ax, Interpreter='latex');

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
  ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
  ax.ZAxis.MinorTick = 'on';
  ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 10.0];
  ax.YAxis.Limits = [0.0, 8.0];
  ax.ZAxis.Limits = [0.0, 18.0];

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

  % Grid lines
  % ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

  % 3D plot view
  view(45, 15.0);

end