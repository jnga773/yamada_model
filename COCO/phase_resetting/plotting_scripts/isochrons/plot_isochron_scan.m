function plot_isochron_scan(isochron_run_in)
  % plot_isochron_scan(run_in)
  %
  % Plots the isochrons from each subdirectory in 'run_in'.

  %-----------------------------%
  %     Read Data: Isochrons    %
  %-----------------------------%
  % Folder name
  dir_data = sprintf('./data/%s/', isochron_run_in);
  % List all directories
  dirs = dir(dir_data);
  % Remove ./ and ../
  dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
  % Sub folder names
  dir_sub = {dirs.name};

  % Empty arrays for theta_old and theta_new
  iso1_data = {};
  iso2_data = {};
  iso3_data = {};
  
  dir_sub_plot = dir_sub;

  % Cycle through each data directory
  for i = 1 : length(dir_sub_plot)
    % Sub folder name
    dir_read = {isochron_run_in, dir_sub_plot{i}};
    fprintf('run_dir:  {%s, %s} \n', dir_read{1}, dir_read{2});

    % Bifurcation data
    bd_read = coco_bd_read(dir_read);

    % Read theta_old and theta_new
    iso1_read = coco_bd_col(bd_read, 'iso1');
    iso2_read = coco_bd_col(bd_read, 'iso2');
    iso3_read = coco_bd_col(bd_read, 'iso3');

    % Append to arrays
    iso1_data{i} = iso1_read;
    iso2_data{i} = iso2_read;
    iso3_data{i} = iso3_read;

  end

  %--------------------------------------%
  %     Read Data: Unperturbed Orbit     %
  %--------------------------------------%
  % Read unperturbed periodic orbit data
  load('./data_mat/initial_PO.mat');

  %-------------------------------------------------------------------------%
  %%                         Plot: Multi Isochrons                         %%
  %-------------------------------------------------------------------------%
  % matplotlib colour order
  colours = colororder();

  fig = figure(1);
  fig.Name = 'Single (Test) Isochron';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 8]; fig.PaperSize = [8, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  hold(ax, 'on');

  % Plot unperturbed periodic orbit
  plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

  % Plot stable manifold of q / x_{+}
  plot3(ax, W_q_stable(:, 1), W_q_stable(:, 2), W_q_stable(:, 3), ...
        Color=colours(1, :), ...
        DisplayName='$W^{s}(p)$');

  % Cycle through data and plot
  for i = 1 : length(dir_sub_plot)

    % Read data
    iso1_plot = iso1_data{i};
    iso2_plot = iso2_data{i};
    iso3_plot = iso3_data{i};

    % Plot single isochron
    plot3(ax, iso1_plot, iso2_plot, iso3_plot, ...
          Color=[0, 0, 0, 0.5], LineStyle='-', ...
          HandleVisibility='off');

  end

  % % Read data
  % i = 1;
  % iso1_plot = iso1_data{1};
  % iso2_plot = iso2_data{1};
  % iso3_plot = iso3_data{1};
  % 
  % % Plot single isochron
  % plot3(ax, iso1_plot, iso2_plot, iso3_plot, ...
  %       Color=[0, 0, 0, 0.5], LineStyle='-', ...
  %       HandleVisibility='off');

  % Legend
  legend(ax, Interpreter='latex');

  % Turn of axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % % X-Axis
  % ax.XAxis.TickValues = 0.0 : 1.0 : 5.0;
  % ax.XAxis.MinorTick = 'on';
  % ax.XAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

  % % Y-Axis
  % ax.YAxis.TickValues = 0.0 : 1.0 : 5.0;
  % ax.YAxis.MinorTick = 'on';
  % ax.YAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

  % % Z-Axis
  % ax.ZAxis.TickDirection = 'in';
  % % ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
  % % ax.ZAxis.MinorTick = 'on';
  % % ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 6.0];
  ax.YAxis.Limits = [0.0, 6.0];
  % ax.ZAxis.Limits = [0.0, ceil(max(xbp_PO(:, 3)))];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$G(t)$';
  ax.YAxis.Label.String = '$Q(t)$';
  ax.ZAxis.Label.String = '$I(t)$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Initial Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % 3D plot view
  view(45, 15.0);


end