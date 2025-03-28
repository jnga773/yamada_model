% Clear plots
% close all;
function plot_bifurcation_diagram(run_names_in)
  % Run names
  H_run = run_names_in.hopf_bifurcations;
  S_run = run_names_in.saddle_nodes;
  T_run = run_names_in.transcritical;
  L_run = run_names_in.approx_homo.continue_homoclinics;
  % L_run = run_names_in.lins_method.continue_homoclinics;
  D_run = run_names_in.limit_cycle.follow_limit_cycle;
  
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read COCO data matrices
  bd_H = coco_bd_read(H_run);
  bd_S = coco_bd_read(S_run);
  bd_T = coco_bd_read(T_run);
  bd_L = coco_bd_read(L_run);
  bd_D = coco_bd_read(D_run);
  
  % Hopf bifurcation line (H)
  A_run8 = coco_bd_col(bd_H, 'A');
  gamma_run8 = coco_bd_col(bd_H, 'gamma');

  % Find the minimum to split into H line and NSA line
  [~, idx] = min(A_run8);
  A_NSA = A_run8(1:idx); gamma_NSA = gamma_run8(1:idx);
  A_H = A_run8(idx+1:end); gamma_H = gamma_run8(idx+1:end);

  % Saddle-Node bifurcation line (A_S)
  A_SN = coco_bd_col(bd_S, 'A');
  gamma_SN = coco_bd_col(bd_S, 'gamma');

  % Transcritical bifurcation line (A_T)
  A_T = coco_bd_col(bd_T, 'A');
  gamma_T = coco_bd_col(bd_T, 'gamma');

  % Approximate homoclinic line
  A_L = coco_bd_col(bd_L, 'A');
  gamma_L = coco_bd_col(bd_L, 'gamma');

  % Approximate double limit cycle line
  A_D = coco_bd_col(bd_D, 'A');
  gamma_D = coco_bd_col(bd_D, 'gamma');

  %-------------------------%
  %     Read Parameters     %
  %-------------------------%
  % Read B and a parameters
  B = coco_bd_val(bd_H, 1, 'B');
  a = coco_bd_val(bd_H, 1, 'a');

  %-----------------------------------------------------------------------%
  %                          Plot: Big Picture                            %
  %-----------------------------------------------------------------------%
  % Plot colours
  colours = colororder();

  %----------------------%
  %     Figure Setup     %
  %----------------------%
  fig = figure(1); clf;
  fig.Name = 'Yamada-Bifurcations';
  fig.Units = 'inches';
  fig.Position = [3, 3, 12, 8];
  fig.PaperSize = [12, 8];

  % Setup axis
  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  % Fontsize
  ax.FontSize = 14;

  % Turn on axis hold
  hold(ax, 'on');

  %-------------------------------%
  %     Plot: MATLAB Plotting     %
  %-------------------------------%
  % Plot saddle-node bifurcations from run4
  plot(ax, A_SN, gamma_SN, LineStyle='-', DisplayName='S');

  % Plot transcritical bifurcations from run5
  plot(ax, A_T, gamma_T, LineStyle='--', DisplayName='T');

  % Plot Hopf birfurcations from run3
  plot(ax, A_H, gamma_H, LineStyle='-', DisplayName='H');

  % Plot Neutral Saddle-Node line run3
  plot(ax, A_NSA, gamma_NSA, LineStyle='--', DisplayName='NS')

  % Plot approximate homoclinic from run8
  plot(ax, A_L, gamma_L, LineStyle='-', DisplayName='L');

  % Plot double limit cycle line from run11
  plot(ax, A_D, gamma_D, LineStyle='-', DisplayName='D');

  % Legend
  legend(ax, 'Interpreter', 'latex');

  % Turn off axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickDirection = 'in';
  ax.XAxis.TickValues = 5.0:1.0:12.0;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 5.0:0.2:12.0;
  
  % Y-Axis
  ax.YAxis.TickValues = 0.0:0.05:0.30;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = 0.00:0.01:0.30;

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$A$';
  ax.YAxis.Label.String = '$\gamma$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  title_str = sprintf('Yamada Model with $\\left( B = %.2f, a = %.2f \\right)$', B, a);
  ax.Title.String = title_str;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [5.0, 11.0];
  ax.YAxis.Limits = [0.0, 0.25];

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  exportgraphics(fig, './images/yamada_bifurcations.pdf', ContentType='vector');

end