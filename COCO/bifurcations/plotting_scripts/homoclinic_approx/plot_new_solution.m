function plot_new_solution(run_old_in, run_new_in, save_figure)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read solution of previous run with largest period, in this case LABEL=142
  % label_max = coco_bd_labs(bd_old, 'EP');
  % label_max = max(label_max);
  [sol_old, ~] = po_read_solution('', run_old_in, max(coco_bd_labs(coco_bd_read(run_old_in), 'EP')));
  [sol_new, ~] = po_read_solution('homo', run_new_in, 2);

  % x solution
  x0_old = sol_old.xbp;
  x0_new = sol_new.xbp;

  %-----------------------------------------------------------------------%
  %                        Plot: New Solution (2D)                        %
  %-----------------------------------------------------------------------%
  fig = figure(9); clf;
  fig.Name = 'New Solution Periodic Orbit (2D)';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [8, 6];

  % Setup axes
  tiles = tiledlayout(1, 3, Padding='compact', TileSpacing='compact');
  ax1 = nexttile; ax2 = nexttile; ax3 = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax1, 'on');
  hold(ax2, 'on');
  hold(ax3, 'on');

  % Plot: run1_in solutions
  leg_lab = sprintf('$\\left( T = %.2f \\right)$', sol_old.T);
  plot(ax1, x0_old(:, 1), x0_old(:, 2), LineWidth=2.0, ...
      DisplayName=leg_lab);
  plot(ax2, x0_old(:, 1), x0_old(:, 3), LineWidth=2.0, ...
      DisplayName=leg_lab);
  plot(ax3, x0_old(:, 2), x0_old(:, 3), LineWidth=2.0, ...
      DisplayName=leg_lab);

  % Plot: run2_in solutions
  leg_lab = sprintf('$\\left( T = %.2f \\right)$', sol_new.T);
  plot(ax1, x0_new(:, 1), x0_new(:, 2), LineStyle='--', Marker='none', ...
      DisplayName=leg_lab);
  plot(ax2, x0_new(:, 1), x0_new(:, 3), LineStyle='--', Marker='none', ...
      DisplayName=leg_lab);
  plot(ax3, x0_new(:, 2), x0_new(:, 3), LineStyle='--', Marker='none', ...
      DisplayName=leg_lab);

  % Legend
  legend(ax1, 'Interpreter', 'latex')
  legend(ax2, 'Interpreter', 'latex')
  legend(ax3, 'Interpreter', 'latex')

  hold(ax1, 'off');
  hold(ax2, 'off');
  hold(ax3, 'off');

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax1.XAxis.Limits = [0.0, 7.0];
  ax1.YAxis.Limits = [0.0, 5.5];

  ax2.XAxis.Limits = [0.0, 7.0];
  ax2.YAxis.Limits = [0.0, 6.5];

  ax3.XAxis.Limits = [0.0, 7.0];
  ax3.YAxis.Limits = [0.0, 6.5];

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
  ax2.Title.String = 'High Period Periodic Orbit (New Solution)';

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
    exportgraphics(fig, './images/new_solution.pdf', ContentType='vector');
  end

end