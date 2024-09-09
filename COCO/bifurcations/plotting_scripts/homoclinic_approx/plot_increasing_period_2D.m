function plot_increasing_period_2D(run_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Grab the label of the run with the highest period
  this_run = run_in;
  bd_new = coco_bd_read(this_run);
  label_max = coco_bd_labs(bd_new, 'EP');
  label_max = max(label_max);

  %-----------------------------------------------------------------------%
  %                    Plot: Increasing Period PO (2D)                    %
  %-----------------------------------------------------------------------%
  fig = figure(4); clf;
  fig.Name = 'Sample-Periodic Orbits (2D)';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [8, 6];

  % Setup subplots
  tiles = tiledlayout(1, 3, Padding='compact', TileSpacing='compact');
  ax1 = nexttile; ax2 = nexttile; ax3 = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  hold(ax1, 'on');
  hold(ax2, 'on');
  hold(ax3, 'on');

  % Loop through indices in [1:1:10, 15:5:100]
  for i = [1:2:label_max-1, label_max]
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
      plot(ax1, xbp_plot(:, 1), xbp_plot(:, 2), Color=colour, LineWidth=lw);
      plot(ax2, xbp_plot(:, 1), xbp_plot(:, 3), Color=colour, LineWidth=lw);
      plot(ax3, xbp_plot(:, 2), xbp_plot(:, 3), Color=colour, LineWidth=lw);
  end

  hold(ax1, 'off');
  hold(ax2, 'off');
  hold(ax3, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax1.XAxis.Label.String = '$G$';
  ax1.YAxis.Label.String = '$Q$';

  ax2.XAxis.Label.String = '$G$';
  ax2.YAxis.Label.String = '$I$';

  ax3.XAxis.Label.String = '$Q$';
  ax3.YAxis.Label.String = '$I$';
  
  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax2.Title.String = 'High Period Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax1, 'on'); grid(ax1, 'on');
  box(ax2, 'on'); grid(ax2, 'on');
  box(ax3, 'on'); grid(ax3, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  % if save_figure == true
  %   exportgraphics(fig, './images/2D_increasing_period.pdf', ContentType='vector');
  % end

end