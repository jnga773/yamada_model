function plot_phase_transition_curve(run_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Bifurcation data
  bd = coco_bd_read(run_in);

  % theta_old data
  theta_old = coco_bd_col(bd, 'theta_old');
  % theta_new data
  theta_new = coco_bd_col(bd, 'theta_new');

  % Read A_perturb value
  A_perturb = coco_bd_val(bd, 1, 'A_perturb');
  % Read theta_perturb and phi_perturb
  theta_perturb = coco_bd_val(bd, 1, 'theta_perturb');
  phi_perturb   = coco_bd_val(bd, 1, 'phi_perturb');

  % Directional vector
  d_vec = [cos(theta_perturb);
           0.0;
           cos(theta_perturb)];

  %--------------%
  %     Plot     %
  %--------------%
   % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(4); clf;
  fig.Name = 'Theta Bifurcations';
  fig.Units = 'inches'; fig.Position = [0, 0, 8, 8]; fig.PaperSize = [8, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  hold(ax, 'on');

  % Plot phase resetting curve
  plot(ax, theta_old, theta_new, LineStyle='-', Color=colours(1, :));

  % Plot diagonal line
  plot(ax, [0, 1], [0, 1], LineStyle='--', Color=[0, 0, 0, 0.75]);

  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  ax.XAxis.TickValues = 0.0 : 0.1 : 1.0;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 0.05 : 0.1 : 1.0;

  ax.YAxis.TickValues = -0.1 : 0.1 : 2.1;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = -0.5 : 0.1 : 2.1;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [-0.025, 1.025];
  ax.YAxis.Limits = [-0.1, 2.1];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$\theta_{\mathrm{old}}$';
  ax.YAxis.Label.String = '$\theta_{\mathrm{new}}$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  title_str = sprintf('Phase Transition Curve (PTC) with $A = %.2f$ and $\\vec{d} = (%.0f, %.0f, %.0f)$', A_perturb, d_vec(1), d_vec(2), d_vec(3));
  ax.Title.String = title_str;

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %----------------------%
  %      Save Figure     %
  %----------------------%
  % Filename
  figname = 'PTC_single';
  exportgraphics(fig, ['./images/', figname, '.pdf'], ContentType='vector');
  
end
