 function plot_initial_periodic_orbits(t_plot, x_plot, run_in, save_figure)
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
  % this_run = run9;
  % label_max = max(coco_bd_labs(coco_bd_read(this_run), 'EP'));
  % % label_max = 1;
  % [sol_plot, ~] = coll_read_solution('po.orb', this_run, label_max);
  % 
  % % Read state and parameters from solution
  % x_soln = sol_plot.xbp;
  % t_soln = sol_plot.tbp;

  %-----------------------------------------------------------------------%
  %                            Plot: Sample PO                            %
  %-----------------------------------------------------------------------%
  fig = figure(10); clf;
  fig.Name = 'Sample-Periodic Orbits';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [6, 8];

  % ax1 = subplot(3, 1, 1); ax2 = subplot(3, 1, 2); ax3 = subplot(3, 1, 3);

  tiles = tiledlayout(1, 3, Padding='compact', TileSpacing='compact');
  ax1 = nexttile; ax2 = nexttile; ax3 = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax1, 'on');
  hold(ax2, 'on');
  hold(ax3, 'on');

  % Plot solutions: To steady state
  plot(ax1, t_plot, x_plot(:, 1), LineStyle='-')
  plot(ax2, t_plot, x_plot(:, 2), LineStyle='-')
  plot(ax3, t_plot, x_plot(:, 3), LineStyle='-')

  % Turn of axis hold
  hold(ax1, 'off');
  hold(ax2, 'off');
  hold(ax3, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax1.XAxis.Label.String = '$t$';
  ax1.YAxis.Label.String = '$G(t)$';

  ax2.XAxis.Label.String = '$t$';
  ax2.YAxis.Label.String = '$Q(t)$';

  ax3.XAxis.Label.String = '$t$';
  ax3.YAxis.Label.String = '$I(t)$';
  
  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax2.Title.String = 'Initial Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax1, 'on'); grid(ax1, 'on');
  box(ax2, 'on'); grid(ax2, 'on');
  box(ax3, 'on'); grid(ax3, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  if save_figure == true
    exportgraphics(fig, './images/initial_periodic_orbits.pdf', ContentType='vector');
  end

end