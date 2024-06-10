function test_compare_homoclinics(run_names_in, save_figure)
  % Run names
  run_approx = run_names_in.approx_homo.continue_homoclinics;

  run_lins   = 'homoclinic_wut';
  run_runs   = {'run1', 'run2', 'run3', 'run4'};
  % run_runs = {'run2'};

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read COCO data matrices
  bd_approx = coco_bd_read(run_approx);
  % Approximate homoclinic data
  A_approx     = coco_bd_col(bd_approx, 'A');
  gamma_approx = coco_bd_col(bd_approx, 'gamma');

  % Lin's method data
  % Cycle through data files
  A_lins = {}; gamma_lins = {};
  for i = 1:length(run_runs)
    % Bifurcation data structure
    homo_run = {run_lins, run_runs{i}};
    bd_lins = coco_bd_read(homo_run);

    % Read data
    A_read     = coco_bd_col(bd_lins, 'A');
    gamma_read = coco_bd_col(bd_lins, 'gamma');

    % Append data
    A_lins{i}     = A_read;
    gamma_lins{i} = gamma_read;
  end

  %-----------------------------------------------------------------------%
  %                          Plot: Big Picture                            %
  %-----------------------------------------------------------------------%
  % Plot colours
  colours = colororder();

  %----------------------%
  %     Figure Setup     %
  %----------------------%
  fig = figure(1); fig.Name = 'Yamada-Bifurcations'; clf;
  fig.Units = 'inches'; fig.Position = [0, 0, 12, 8]; fig.PaperSize = [12, 8];
  ax = gca();
  ax.FontSize = 14.0; 

  %--------------%
  %     Plot     %
  %--------------%
  % Turn on axis hold
  hold(ax, 'on');

  % Plot approximate homoclinic from run8
  plot(ax, A_approx, gamma_approx, LineStyle='-', DisplayName='Homoclinic (Approximate)');

  % Colour for second plot
  C1 = [255, 127,  14] / 255;
  for i = 1:length(run_runs)
    if i == 1
      plot(ax, A_lins{i}, gamma_lins{i}, LineStyle='--', Color=C1, ...
           DisplayName="Homoclinc (Lin's)");
    else
      plot(ax, A_lins{i}, gamma_lins{i}, LineStyle='--', Color=C1, ...
           HandleVisibility='off');
    end
  end

  % % Plot MX points
  % A_MX     = [6.0763, 6.2376, 6.2571, 6.6042, 6.6411, 6.7329];
  % gamma_MX = [1.9371e-1, 1.4565e-1, 1.4163e-1, 8.8593e-2, 8.3900e-2, 7.1591e-2];
  % plot(ax, A_MX, gamma_MX, LineStyle='none', Color='r', Marker='x', ...
  %      MarkerSize=15, DisplayName='MX');

  % Legend
  legend(ax, 'Interpreter', 'latex');

  % Turn off axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickValues = 6.0:0.1:6.9;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 6.05:0.1:6.85;

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
  title_str = sprintf('Yamada Model (Comparing Homoclinic Methods)');
  ax.Title.String = title_str;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [6.0, 6.9];
  ax.YAxis.Limits = [0.0, 0.25];

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  if save_figure == true
    exportgraphics(fig, './images/compared_approx_and_lins.pdf', ContentType='vector');
  end

end