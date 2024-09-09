function plot_extended_time_profile(run_in, label_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  sol_plot = po_read_solution('homo', run_in, label_in);

  % Evaluate vector field at basepoints
  f = yamada(sol_plot.xbp', repmat(sol_plot.p, [1 size(sol_plot.xbp, 1)]));

  % Find basepoint closest to equilibrium
  [~, idx] = min(sqrt(sum(f.*f, 1))); 

  % Reconstruct some time bidniz
  t_plot = [sol_plot.tbp(idx:end) - sol_plot.tbp(idx); ...
            sol_plot.tbp(1:idx) + (sol_plot.T - sol_plot.tbp(idx))];

  x_plot = sol_plot.xbp([idx:end, 1:idx], 3);

  %-----------------------------------------------------------------------%
  %                      Plot: Extended Time Profile                      %
  %-----------------------------------------------------------------------%
  fig = figure(5); clf;
  fig.Name = 'Extended Time Profile';
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
  % plot(tbp, xbp, 'LineStyle', '-', 'LineWidth', 2, ...
  %      'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  plot(ax, t_plot, x_plot, LineStyle='-', LineWidth=2, Color='black', ...
      Marker='.', MarkerSize=12)
  
  % Turn off axis hold
  hold(ax, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$t$';
  ax.YAxis.Label.String = '$x_{3}(t)$';
  
  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Extended Time Profile';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  % if save_figure == true
  %   exportgraphics(fig, './images/extended_time_profile.pdf', ContentType='vector');
  % end

end
