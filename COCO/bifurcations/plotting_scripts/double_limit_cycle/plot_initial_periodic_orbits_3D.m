 function plot_initial_periodic_orbits_3D(t_plot, x_plot, run_in, save_figure)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Long-time data
  % t_plot = t_long;
  % x_plot = x_long;

  % Steady state data
  % t_plot = t_PO;
  % x_plot = x_PO;

  % Read solution from COCO run
  this_run = run_in;
  label_max = max(coco_bd_labs(coco_bd_read(this_run), 'EP'));
  % label_max = 1;
  % label_max = 33;
  [sol_plot, ~] = coll_read_solution('po.orb', this_run, label_max);

  % Read state and parameters from solution
  x_soln = sol_plot.xbp;

  %-----------------------------------------------------------------------%
  %                         Plot: Initial PO (3D)                         %
  %-----------------------------------------------------------------------%
  % Plot colours
  colours = colororder();

  %----------------------%
  %     Figure Setup     %
  %----------------------%
  fig = figure(12); clf;
  fig.Name = 'Sample-Periodic Orbits (2D)';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [8, 6];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plot solutions: To steady state
  plot3(ax, x_plot(:, 1), x_plot(:, 2), x_plot(:, 3), LineStyle='-')

  % Plot solutions: From COCO
  plot3(ax, x_soln(:, 1), x_soln(:, 2), x_soln(:, 3), LineStyle='--')

  % Turn of axis hold
  hold(ax, 'off');

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
  box(ax, 'on'); grid(ax, 'on');
  
  % 3D view
  view(45, 15.0);

  %---------------------%
  %     Save Figure     %
  %---------------------%
  if save_figure == true
    % exportgraphics(fig, './images/initial_periodic_orbits_3D.png', Resolution=800);
    exportgraphics(fig, './images/initial_periodic_orbits_3D.pdf', ContentType='vector');
  end
end