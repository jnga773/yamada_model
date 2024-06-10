function plot_A_vs_T(run_in, save_figure)
  %-------------------%
  %     Read Data     %
  %-------------------%
  A_data = coco_bd_col(coco_bd_read(run_in), 'A')';
  T_data = coco_bd_col(coco_bd_read(run_in), 'po.period')';

  %-----------------------------------------------------------------------%
  %                        Plot A vs. T from run6                         %
  %-----------------------------------------------------------------------%
  fig = figure(3); clf;
  fig.Name = 'A vs. T';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [8, 6];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axis
  hold(ax, 'on');

  % Plot
  plot(A_data, T_data, LineStyle='-', LineWidth=1.0)
  
  % Turn off axis hold
  hold(ax, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$A$';
  ax.YAxis.Label.String = '$T$';
  
  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Pump Current $(A)$ vs Period $(T)$';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  if save_figure == true
    exportgraphics(fig, './images/A_vs_T.pdf', ContentType='vector');
  end

end
