function plot_single_DTC(run_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  bd_read = coco_bd_read(run_in);

  % theta_old
  theta_old = coco_bd_col(bd_read, 'theta_old');
  % theta_new
  theta_new = coco_bd_col(bd_read, 'theta_new');
  % theta_perturb
  theta_perturb = coco_bd_col(bd_read, 'theta_perturb');
  % A_perturb
  A_perturb = coco_bd_col(bd_read, 'A_perturb');

  %-------------------%
  %     Plot Data     %
  %-------------------%
  colours = colororder();

  fig = figure(1); clf;
  ax = gca();

  width = 6.0;
  height = 8.0;

  set_figure_dimensions(width, height);
  

  % Plot
  hold(ax, 'on');

  % Fundamental domain
  patch([-3, 3, 3, -3], [0, 0, 1, 1], colours(3, :), ...
        FaceAlpha=0.2, EdgeColor='none');
  

  % DTC
  plot(ax, theta_perturb / pi, theta_new, LineStyle='-', Color=colours(1, :));

  % daspect(ax, [1, 1, 1]);

  % Limits
  xlim(ax, [0.0, 0.5]);
  ylim(ax, [-0.1, 1.1]);

  % Labels
  xlabel(ax, '$\varphi_{\mathrm{d}} / \pi$');
  ylabel(ax, '$\vartheta_{\mathrm{n}}$');
  title_str = sprintf('$\\vartheta_{\\mathrm{o}} = %.4f, A_{\\mathrm{p}} = %.4f$', theta_old(1), A_perturb(1));
  title(ax, title_str);

  % Figure stuff
  box(ax, 'on');

end