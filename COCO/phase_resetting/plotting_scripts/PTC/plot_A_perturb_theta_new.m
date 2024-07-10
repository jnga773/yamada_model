function plot_A_perturb_theta_new(run_in)
  % plot_A_perturb_theta_new(run_in)
  %
  % Plots the perturbation amplitude against the resulting \theta_new.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Bifurcation data
  bd = coco_bd_read(run_in);

  % Read perturbation amplitudes
  A_perturb = coco_bd_col(bd, 'A_perturb');
  % Read theta_new
  theta_new = coco_bd_col(bd, 'theta_new');

  % Read perturbation vector angles
  theta_perturb = coco_bd_val(bd, 1, 'theta_perturb');
  phi_perturb   = coco_bd_val(bd, 1, 'phi_perturb');

  % Perturbation directional vector
  d_vec = [cos(theta_perturb) * sin(phi_perturb);
           0.0;
           sin(theta_perturb)];

  % Get all SP labels
  SP_labs = coco_bd_labs(bd, 'SP');

  % Read perturbation and theta_new values for each SP point
  SP_A_perturb = zeros(length(SP_labs), 1);
  SP_theta_new = zeros(length(SP_labs), 1);

  for i = 1 : length(SP_labs)
    SP_A_perturb(i) = coco_bd_val(bd, SP_labs(i), 'A_perturb');
    SP_theta_new(i) = coco_bd_val(bd, SP_labs(i), 'theta_new');
  end

  %-----------------------------------------------------------------------%
  %                   Plot A_perturb against theta_new                    %
  %-----------------------------------------------------------------------%
  % Plotting colours
  colours = colororder();

  fig = figure(3); clf;
  fig.Name = 'PTC Scans';
  fig.Units = 'inches'; fig.Position = [0, 0, 8, 8]; fig.PaperSize = [8, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;
  ax.FontSize = 11;
  
  %--------------%
  %     Plot     %
  %--------------%
  hold(ax, 'on');
  plot(ax, A_perturb, theta_new);
  plot(ax, SP_A_perturb, SP_theta_new, LineStyle='none', Marker='o', MarkerSize=14, ...
       MarkerEdgecolor=colours(2, :), MarkerFaceColor=colours(2, :))
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % ax.XAxis.TickValues = 0.0 : 0.1 : 1.0;
  % ax.XAxis.MinorTick = 'on';
  % ax.XAxis.MinorTickValues = 0.05 : 0.1 : 1.0;

  ax.YAxis.TickValues = -0.1 : 0.1 : 1.1;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = -0.05 : 0.1 : 1.1;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [-0.05, max(A_perturb) + 0.05];
  ax.YAxis.Limits = [-0.01, 1.025];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$A_{\mathrm{perturb}}$';
  ax.YAxis.Label.String = '$\theta_{\mathrm{new}}$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  title_str = sprintf('Phase Transition Curve (PTC) with $\\vec{d} = (%.0f, %.0f, %.0f)$', d_vec(1), d_vec(2), d_vec(3));
  ax.Title.String = title_str;

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % %---------------------%
  % %     Save Figure     %
  % %---------------------%
  % % Filename
  % exportgraphics(fig, './images/A_perturb_theta_new.pdf', ContentType='vector');

end