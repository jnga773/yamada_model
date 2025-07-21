function plot_perturbed_states(run_in, label_in)
  % plot_perturbed_states(run_in, label_in)
  %
  % Plots the perturbed orbit segments in time

  %--------------------------------------%
  %     Read Data: Perturbed Segment     %
  %--------------------------------------%
  % Read segment data
  sol4 = coll_read_solution('seg4', run_in, label_in);

  % State space data
  x4   = sol4.xbp;

  % Time data
  t4   = sol4.tbp;

  % Parameters
  p = sol4.p;
  T = p(5);
  k = p(6);
  theta_old = p(7);

  % Renormalise time data
  t4 = t4 * k;

  % Add \theta_old to time array for the point where the perturbation
  % is applied
  if theta_old < 1.0
    t4 = t4 + theta_old;
  else
    t4 = t4 + mod(theta_old, 1.0);
  end

  %--------------------------------------%
  %     Read Data: Unperturbed Orbit     %
  %--------------------------------------%
  % Load in initial periodic orbit data
  load('./data_mat/initial_PO.mat', 'xbp_PO', 'tbp_PO');

  % append x_PO and t_PO 25 times
  x_PO_plot = [];
  t_PO_plot = [];
  
  tbp_PO = tbp_PO / max(tbp_PO);
  t_max = 0.0;

  for i = 1 : k
    % Append x_PO
    x_PO_plot = [x_PO_plot; xbp_PO(2:end, :)];

    % Append t_PO
    t_PO_plot = [t_PO_plot; t_max + tbp_PO(2:end)];

    t_max = t_PO_plot(end);
  end

  % Normalise by k
  % t_PO_plot = t_PO_plot / k;s
      
  %-------------------%
  %     Plot Data     %
  %-------------------%
  fig = figure(1); clf;
  fig.Name = 'Initial Periodic Orbits';
  fig.Units = 'inches'; fig.Position = [3, 3, 14, 10]; fig.PaperSize = [14, 10];

  tiles = tiledlayout(3, 1, Padding='compact', TileSpacing='compact');
  ax1 = nexttile;
  ax2 = nexttile;
  ax3 = nexttile;

  ax = [ax1; ax2; ax3];

  %--------------%
  %     Plot     %
  %--------------%
  for i = 1 : 3
    hold(ax(i), 'on');
    
    % Plot
    plot(ax(i), t4, x4(:, i), LineStyle='-', LineWidth=2.5)
    plot(ax(i), t_PO_plot, x_PO_plot(:, i), LineStyle=':')

    hold(ax(i), 'off');
  end

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax3.XAxis.Label.String = '$T$';

  ax1.YAxis.Label.String = '$G(t)$';
  ax2.YAxis.Label.String = '$Q(t)$';
  ax3.YAxis.Label.String = '$I(t)$';

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  for i = 1 : 3
    ax(i).XAxis.Limits = [0.0, k];
  end

  %--------------------%
  %     Axis Title     %
  %--------------------%
  title_str = sprintf('COCO Solution (run: $\\verb!%s!$, label: %d)', run_in, label_in);
  ax1.Title.String = title_str;

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

end
