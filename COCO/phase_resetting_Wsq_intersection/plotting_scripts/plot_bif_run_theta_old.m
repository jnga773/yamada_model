function plot_bif_run_theta_old(run_in)
  % plot_bif_run_theta_old(run_in)
  %
  % Plots the bifurcation diagram of theta_old vs A_perturb.
  %
  % Parameters
  % ----------
  % run_in : str
  %    Run string identifier.
  
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read bifurcation data
  bd_read = coco_bd_read(run_in);

  % Read theta_old values
  theta_old     = coco_bd_col(bd_read, 'theta_gamma');
  % Read A_perturb data
  A_perturb     = coco_bd_col(bd_read, 'A_perturb');
  % Read theta_perturb data
  theta_perturb = coco_bd_col(bd_read, 'theta_perturb');

  % Mod theta_old by 1
  % theta_old = mod(theta_old, 1);

  % % Get unique indices
  % [~, unique_idx] = unique(theta_old);
  % theta_old     = theta_old(unique_idx);
  % A_perturb     = A_perturb(unique_idx);
  % theta_perturb = theta_perturb(unique_idx);

  % Sort by theta_old
  [~, sort_idx] = sort(theta_perturb);
  theta_old     = theta_old(sort_idx);
  A_perturb     = A_perturb(sort_idx);
  theta_perturb = theta_perturb(sort_idx);

  %-------------------%
  %     Plot Data     %
  %-------------------%
  colours = colororder();
  
  % Setup figure
  fig = figure(3); clf;
  fig.Name = 'Periodic Orbit Phase Portrait (3D)';
  ax = gca();

  % Axis dimensions
  width = 7.5;
  height = 4.0;

  % Set figure size
  set_figure_dimensions(width, height, scale=4);

  % Set axis linewidth
  ax.LineWidth = 0.8;

  %-------------------%
  %     Hold Axis     %
  %-------------------%
  hold(ax, 'on');


  %------------------------------------------%
  %     Plot: theta_perturb vs theta_old     %
  %------------------------------------------%
  % Plot: theta_old vs theta_perturb
  plot(ax, theta_perturb, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
  plot(ax, theta_perturb, theta_old-1, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
  plot(ax, theta_perturb+2, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);

  %-------------------%
  %     Hold Axis     %
  %-------------------%
  hold(ax, 'off');


  %---------------------%
  %     Axis Limits     %
  %---------------------%
  % X-Axis: theta_perturb
  ax.XAxis.Limits = [-0.125, 1.0];

  % Y-Axis 1: theta_old
  ax.YAxis.Limits = [0.0, 1.0];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  % X-Axis: theta_perturb
  ax.XAxis.Label.String = '$\varphi_{\mathrm{p}}$';

  % Y-Axis 1: theta_perturb
  ax.YAxis.Label.String = '$\vartheta_{\mathrm{o}}$';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  % grid(ax, 'on');

end