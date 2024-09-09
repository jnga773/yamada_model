function plot_delta_vs_T(run_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Grab the label of the run with the highest period
  this_run = run_in;
  bd_new = coco_bd_read(this_run);
  label_max = coco_bd_labs(bd_new, 'EP');
  label_max = max(label_max);

  % Read period from run6
  % T_data = coco_bd_col(bd6, 'po.period')';
  T_data = [];

  addpath('./functions/field/hardcoded/');

  % Cycle through each run in bd6, read the x0 data, and calculate the
  % minimum part
  delta_data = [];
  for i = 1:label_max
      % Read solution% Read solution with maximum period
      [sol_plot, ~] = coll_read_solution('po.orb', run_in, i);
      
      % Read state and parameters from solution
      x_plot = sol_plot.xbp;
      p_plot = sol_plot.p;

      T = sol_plot.T;
      
      % Evaluate vector field at basepoints
      f = yamada(x_plot', repmat(p_plot, [1, size(x_plot, 1)]));

      % Extract the discretisation points corresponding to the minimum value of
      % the norm of the vector field along the longest-period periodic orbit.
      % Find basepoint closest to equilibrium
      f_norm = sqrt(sum(f .* f, 1));
      [mininima, idx] = min(f_norm);

      % Append data
      T_data = [T_data, T];
      delta_data = [delta_data, mininima];
  end

  %-----------------------------------------------------------------------%
  %                       Plot delta vs. T from run6                      %
  %-----------------------------------------------------------------------%
  fig = figure(2); clf;
  fig.Name = 'Period vs. Min Difference from Equilibrium Point';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 6];
  fig.PaperSize = [8, 6];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  % Set log axis
  ax.YAxis.Scale = 'log';

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axis
  hold(ax, 'on');

  % Plot
  plot(T_data, delta_data, LineStyle='-', LineWidth=1.0)
  
  % Turn off axis hold
  hold(ax, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$T$';
  ax.YAxis.Label.String = '$|x_{\mathrm{sol}} - x_{0}|$';
  
  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Period vs. Min Difference from Equilibrium Point';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  % if save_figure == true
  %   exportgraphics(fig, './images/delta_vs_T.pdf', ContentType='vector');
  % end

end