function plot_two_axis_figure(run1_in, run2_in, runq_in)
  % plot_two_axis_figure(run1_in, run2_in)
  %
  % Plots the bifurcation diagram of theta_old vs A_perturb for two runs.
  %
  % Parameters
  % ----------
  % run1_in : str
  %    First run string identifier.
  % run2_in : str
  %    Second run string identifier.
  % runq_in : str
  %    Third run string identifier.

  %-------------------------%
  %     Read Data: Easy     %
  %-------------------------%
  % Read bifurcation data
  bd_readq = coco_bd_read(runq_in);

  % Get labels
  qlabs = coco_bd_labs(bd_readq, 'q_PT');

  % Read theta_old
  theta_oldq = coco_bd_val(bd_readq, 'q_PT', 'theta_gamma')';

  % Read xpos solution
  [solq, ~] = ep_read_solution('xpos', runq_in, qlabs(1));
  xpos = solq.x';

  % Read seg3 solution
  [sol2_1, ~] = coll_read_solution('seg2', runq_in, qlabs(1));
  [sol2_2, ~] = coll_read_solution('seg2', runq_in, qlabs(2));

  % Read state space solutions
  xbp1 = sol2_1.xbp(1, :);
  xbp2 = sol2_2.xbp(1, :);

  % Calculate distances
  vec1 = xpos - xbp1;
  vec2 = xpos - xbp2;

  % Displacement amplitude
  A1 = norm(vec1);
  A2 = norm(vec2);

  % Angle of displacement vector
  T1 = mod(atan2(vec1(3), vec1(1)), 2*pi);
  T2 = mod(atan2(vec2(3), vec2(1)), 2*pi);
  
  % Normalise
  T1 = T1 / (2 * pi);
  T2 = T2 / (2 * pi);

  % Make into arrays
  A_perturbq    = [A1; A2];
  theta_perturbq = [T1; T2];

  %-------------------------%
  %     Read Data: Easy     %
  %-------------------------%
  % Read bifurcation data for both runs
  bd_read1 = coco_bd_read(run1_in);
  bd_read2 = coco_bd_read(run2_in);

  % Read theta_old values
  theta_old1     = coco_bd_col(bd_read1, 'theta_gamma')';
  theta_old2     = coco_bd_col(bd_read2, 'theta_gamma')';
  
  % Read A_perturb data
  A_perturb1     = coco_bd_col(bd_read1, 'A_perturb')';
  A_perturb2     = coco_bd_col(bd_read2, 'A_perturb')';
  
  % Read theta_perturb data
  theta_perturb1 = coco_bd_col(bd_read1, 'theta_perturb')';
  theta_perturb2 = coco_bd_col(bd_read2, 'theta_perturb')';

  %----------------------%
  %     Combine Data     %
  %----------------------%
  % Combine data arrays
  theta_old     = [theta_old1; theta_old2; theta_oldq];
  A_perturb     = [A_perturb1; A_perturb2; A_perturbq];
  theta_perturb = [theta_perturb1; theta_perturb2; theta_perturbq];

  % Mod the data
  theta_perturb = mod(theta_perturb, 1.0);

  %-----------------------------------------%
  %     Sort Data: Increasing theta_old     %
  %-----------------------------------------%
  % Sort
  [~, sort_idx] = sort(theta_old);

  % Sort by theta_old
  theta_old     = theta_old(sort_idx);
  theta_perturb = theta_perturb(sort_idx);
  A_perturb     = A_perturb(sort_idx);
  
  % Sort
  [~, sort_idx] = sort(theta_perturb);

  % Sort by theta_old
  theta_old     = theta_old(sort_idx);
  theta_perturb = theta_perturb(sort_idx);
  A_perturb     = A_perturb(sort_idx);

  % %------------------------------------------%
  % %     Shift Data: By max theta_perturb     %
  % %------------------------------------------%
  % % Find max point of theta_old
  % [max_val, max_idx] = max(theta_perturb);
  % 
  % % Shift data around
  % theta_perturb = [theta_perturb(1:max_idx)-2; theta_perturb(max_idx+1:end)];
  % theta_old     = [theta_old(1:max_idx); theta_old(max_idx+1:end)];
  % A_perturb     = [A_perturb(1:max_idx); A_perturb(max_idx+1:end)];

  %---------------------%
  %     Extend Data     %
  %---------------------%
  % theta_perturb = [theta_perturb-2; theta_perturb; theta_perturb+2];
  % theta_old     = [theta_old-1; theta_old; theta_old+1];
  % A_perturb     = [A_perturb; A_perturb; A_perturb];

  %-------------------------------------------------------------------------%
  %                         Plot Data: 2D Comparison                        %
  %-------------------------------------------------------------------------%
  colours = colororder();

  % Setup figure
  fig = figure(2); clf;
  fig.Name = 'Intersection with Ws(q) and Perturbed Gamma';
  ax = gca();

  % Axis dimensions
  width = 7.5;
  height = 5.0;

  % Set figure size
  set_figure_dimensions(width, height, scale=2);

  % Set axis linewidth
  ax.LineWidth = 0.8;

  %-------------------%
  %     Hold Axis     %
  %-------------------%
  hold(ax, 'on');

  %--------------------------------------------%
  %     Plot: Highlight theta_perturb Area     %
  %--------------------------------------------%
  % Highlight area over 0.0 <= theta_perturb <= 0.5 pi
  patch(ax(1), [0, 0, 0.25, 0.25], [0, 1, 1, 0], colours(3, :), FaceAlpha=0.2, ...
        EdgeColor='none', HandleVisibility='off')

  %------------------------------------------%
  %     Plot: theta_perturb vs theta_old     %
  %------------------------------------------%
  % Left axis
  yyaxis(ax, 'left');

  % Plot: theta_old vs theta_perturb
  plot(ax(1), theta_perturb, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
  plot(ax(1), theta_perturb, theta_old-1, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
  plot(ax(1), theta_perturb+2, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);

  %------------------------------------------%
  %     Plot: theta_perturb vs A_perturb     %
  %------------------------------------------%
  % Right axis
  yyaxis(ax, 'right');

  % Plot: theta_old vs theta_perturb
  plot(ax, theta_perturb, A_perturb, Color=colours(2, :), LineStyle='-', LineWidth=2.5);
  plot(ax, theta_perturb+2, A_perturb, Color=colours(2, :), LineStyle='-', LineWidth=2.5);

  %-------------------%
  %     Hold Axis     %
  %-------------------%
  hold(ax, 'off');

  %----------------------------%
  %     Axis Ticks: Axis 1     %
  %----------------------------%
  % X-Axis: theta_perturb
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.TickValues = -0.25 : 0.25 : 1.5;
  ax.XAxis.MinorTickValues = -0.25 : 0.125 : 1.5;

  % Y-Axis 1: theta_old
  ax.YAxis(1).MinorTick = 'on';
  ax.YAxis(1).TickValues = 0.0 : 0.25 : 1.0;
  ax.YAxis(1).MinorTickValues = 0.0 : 0.125 : 1.0;

  % Y-Axis 2: A_perturb
  ax.YAxis(2).MinorTick = 'on';
  ax.YAxis(2).TickValues = 0.0 : 3.0 : 15.0;
  ax.YAxis(2).MinorTickValues = 0.0 : 1.0 : 15.0;

  %---------------------%
  %     Tick Labels     %
  %---------------------%
  % ax.XAxis.TickLabels = {};
  % ax.YAxis(1).TickLabels = {};
  % ax.YAxis(2).TickLabels = {};

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  % X-Axis: theta_perturb
  ax.XAxis.Limits = [-0.125, 1.0];

  % Y-Axis 1: theta_old
  ax.YAxis(1).Limits = [0.0, 1.0];

  % Y-Axis 2: A_perturb
  ax.YAxis(2).Limits = [0.0, 15.0];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  % X-Axis: theta_perturb
  ax.XAxis.Label.String = '$\varphi_{\mathrm{p}}$';

  % Y-Axis 1: theta_perturb
  ax.YAxis(1).Label.String = '$\vartheta_{\mathrm{o}}$';

  % Y-Axis 2: A_perturb
  ax.YAxis(2).Label.String = '$A$';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  % grid(ax, 'on');

end