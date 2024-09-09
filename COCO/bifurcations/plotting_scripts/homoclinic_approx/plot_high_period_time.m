function plot_high_period_time(run_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Grab the label of the run with the highest period
  this_run = run_in;
  bd_new = coco_bd_read(this_run);
  label_max = coco_bd_labs(bd_new, 'EP');
  label_max = max(label_max);

  % Read solution with maximum period
  [sol_plot, ~] = coll_read_solution('po.orb', run_in, label_max);

  % Read state and parameters from solution
  x_plot = sol_plot.xbp;
  p_plot = sol_plot.p;
  t_plot = sol_plot.tbp;

  % Evaluate vector field at basepoints
  f = yamada(x_plot', repmat(p_plot, [1, size(x_plot, 1)]));

  % Extract the discretisation points corresponding to the minimum value of
  % the norm of the vector field along the longest-period periodic orbit.
  % Find basepoint closest to equilibrium
  f_norm = sqrt(sum(f .* f, 1)); f_norm = f_norm';
  [~, idx] = min(f_norm);

  %-----------------------------------------------------------------------%
  %                         Plot: High Period PO                          %
  %-----------------------------------------------------------------------%
  fig = figure(3); clf;
  fig.Name = 'High-Period Periodic Orbit (Time)';
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
  % plot(ax, t_plot, x_plot, LineStyle='-', LineWidth=2, Color='black', ...
  %      Marker='.', MarkerSize=12)
  plot(ax, t_plot, x_plot(:, 1), LineStyle='-', DisplayName='$G(t)$')
  plot(ax, t_plot, x_plot(:, 2), LineStyle='--', DisplayName='$Q(t)$')
  plot(ax, t_plot, x_plot(:, 3), LineStyle='-.', DisplayName='$I(t)$')

  % Plot straight line at t(idx), corresponding to point of minimum norm of f
  xline(t_plot(idx), LineStyle=':', Color='Black', LineWidth=1.5, ...
        DisplayName=sprintf('idx = %d', idx))

  % Legend
  legend(ax, 'Interpreter', 'latex')

  % Turn off axis hold
  hold(ax, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$t$';
  ax.YAxis.Label.String = '$\vec{x}(t)$';
  
  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Largest Period Orbit (Time)';

  % Tick params
  ax.XAxis.TickDirection = 'in'; ax.YAxis.TickDirection = 'in';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  % if save_figure == true
  %   exportgraphics(fig, './images/highest_periodic_orbit_time.pdf', ContentType='vector');
  % end

end
