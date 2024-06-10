function plot_phase_reset_phase_space(run_in, label_in, fig_num_in, save_figure, plot_3d)
  % plot_phase_reset_PO(run_in, label_in)
  %
  % Plots the phase resetting periodic orbit from all four segments.

  %-------------------%
  %     Plot Data     %
  %-------------------%
  fig = figure(fig_num_in); clf;
  fig.Name = 'Initial Periodic Orbits';
  fig.Units = 'inches'; fig.Position = [3, 3, 12, 8]; fig.PaperSize = [12, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;
  ax.FontSize = 11;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  if plot_3d == true
    plot_phase_reset_phase_space_template(ax, run_in, label_in)
  else
    plot_2D_phase_reset_phase_space_template(ax, run_in, label_in)
  end

  % Legend
  legend(ax, 'Interpreter', 'latex');

  % Hold axes
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  if plot_3d == true
    % X-Axis
    ax.XAxis.TickDirection = 'in';
    % ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
    % ax.XAxis.MinorTick = 'on';
    % ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;
    
    % Y-Axis
    ax.YAxis.TickDirection = 'in';
    % ax.YAxis.TickValues = 0.0 : 2.0 : 8.0;
    % ax.YAxis.MinorTick = 'on';
    % ax.YAxis.MinorTickValues = 1.0 : 2.0 : 8.0;

    % Z-Axis
    ax.ZAxis.TickDirection = 'in';
    % ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
    % ax.ZAxis.MinorTick = 'on';
    % ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;
  
  else
    % X-Axis
    ax.XAxis.TickDirection = 'in';
    % ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
    % ax.XAxis.MinorTick = 'on';
    % ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;
    
    % Y-Axis
    ax.YAxis.TickDirection = 'in';
    % ax.YAxis.TickValues = 0.0 : 2.0 : 18.0;
    % ax.YAxis.MinorTick = 'on';
    % ax.YAxis.MinorTickValues = 1.0 : 2.0 : 18.0;
  end

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  if plot_3d == true
    ax.XAxis.Limits = [0.0, 10.0];
    ax.YAxis.Limits = [0.0, 8.0];
    ax.ZAxis.Limits = [0.0, 21.0];
  else
    ax.XAxis.Limits = [0.0, 10.0];
    ax.YAxis.Limits = [0.0, 21.0];
  end

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  if plot_3d == true
    ax.XAxis.Label.String = '$G(t)$';
    ax.YAxis.Label.String = '$Q(t)$';
    ax.ZAxis.Label.String = '$I(t)$';
  else
    ax.XAxis.Label.String = '$G(t)$';
    ax.YAxis.Label.String = '$I(t)$';
  end

  %--------------------%
  %     Axis Title     %
  %--------------------%
  title_str = sprintf('COCO Solution (run: $\\verb!%s!$, label: %d)', run_in, label_in);
  % ax.Title.String = title_str;

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % 3D plot view
  if plot_3d == true
    view(45, 15.0);
  end

  %----------------------%
  %      Save Figure     %
  %----------------------%
  if save_figure == true
    % Filename
    figname = sprintf('phase_reset_curve_%s_%d', run_in(1:5), label_in);
    % exportgraphics(fig, ['./images/', figname, '.png'], Resolution=800);
    exportgraphics(fig, ['./images/', figname, '.pdf'], ContentType='vector');
  end

end