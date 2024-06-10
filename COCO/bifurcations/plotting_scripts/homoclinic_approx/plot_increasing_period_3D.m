function plot_increasing_period_3D(run_in, save_figure)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Grab the label of the run with the highest period
  this_run = run_in;
  bd_new = coco_bd_read(this_run);
  label_max = coco_bd_labs(bd_new, 'EP');
  label_max = max(label_max);
  
  %-----------------------------------------------------------------------%
  %                    Plot: Increasing Period PO (3D)                    %
  %-----------------------------------------------------------------------%
  fig = figure(7); clf;
  fig.Name = 'Sample-Periodic Orbits (3D)';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [8, 6];

  % Setup axes
  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Loop through indices in [1:1:10, 15:5:100]
  for i = [1:label_max]
      % fprintf('i = %d\n', i);

      % Read solution
      [sol_plot, ~] = coll_read_solution('po.orb', run_in, i);

      % x solution
      xbp_plot = sol_plot.xbp;

      % Plot colour
      if (i == label_max)
          colour = 'Red';
          lw = 2.0;
      else
          colour = 'Black';
          lw = 0.5;
      end

      % Plot
      plot3(xbp_plot(:, 1), xbp_plot(:, 2), xbp_plot(:, 3), Color=colour, LineWidth=lw);
  end

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
  ax.Title.String = 'High Period Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  if save_figure == true
    exportgraphics(fig, './images/3D_increasing_period.pdf', ContentType='vector');
  end

end