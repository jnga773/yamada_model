clear all; close all; clc;

%---------------------%
%%     Read Data     %%
%---------------------%
% Check to see if 'data_PTC.mat' exists. If not, read data and create
% it
filename_data = './data_PTC.mat';

read_write_data(filename_data);

% Load data
load(filename_data);

%-------------------------------------------------------------------------%
%%      Nonlinear Photonics 2024 Abstract - Phase Transition Curve       %%
%-------------------------------------------------------------------------%
% Plotting colours
colours = colororder();

fig = figure(5); clf;
fig.Name = 'PTC Scans';
fig.Units = 'inches'; fig.Position = [0, 0, 8, 7]; fig.PaperSize = [8, 7];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 18;

%--------------%
%     Plot     %
%--------------%
hold(ax, 'on');

% Cycle through and plot
for i = 1 : length(dir_sub_plot)
  % Theta_old plot
  theta_old = theta_old_data{i};
  theta_new = theta_new_data{i};

  % Get any theta_new values < 0
  idx_neg = (theta_new < 0);
  % Get any theta_new values >= 0 and <= 1
  idx_mid = (theta_new >= 0 & theta_new <= 1);
  % Get any theta_new values > 1
  idx_gt1 = (theta_new > 1 & theta_new <= 2);
  % Get any theta_new values > 1
  idx_gt2 = (theta_new > 2);

  % Legend label
  leg_lab = sprintf('$A_{\\mathrm{p}} = %0.2f$', A_perturb(i));

  % Set colours and line styles
  if i <= 10
    plot_colour = colours(i, :);
    plot_ls     = '-';
  elseif i > 10
    plot_colour = colours(i-10, :);
    plot_ls     = ':';
  end

  %-----------------------------%
  %     Plot: Straight Data     %
  %-----------------------------%
  % Plot: 0 <= \theta_new <= 1
  plot(ax, theta_old(idx_mid), theta_new(idx_mid), ...
       LineStyle=plot_ls, Color=plot_colour, LineWidth=4.0, ...
       DisplayName=leg_lab);
  % Plot: \theta_new < 0
  plot(ax, theta_old(idx_neg), theta_new(idx_neg), ...
       LineStyle='--', Color=plot_colour, LineWidth=4.0, ...
       HandleVisibility='off');
  % Plot: 1 < \theta_new <= 2
  plot(ax, theta_old(idx_gt1), theta_new(idx_gt1), ...
       LineStyle='--', Color=plot_colour, LineWidth=4.0, ...
       HandleVisibility='off');

  %----------------------------%
  %     Plot: Shifted Data     %
  %----------------------------%
  % Plot: \theta_new < 0
  plot(ax, theta_old(idx_neg), theta_new(idx_neg)+1, ...
        LineStyle=plot_ls, Color=plot_colour, LineWidth=4.0, ...
        HandleVisibility='off');
  % Plot: 1 < \theta_new <= 2
  plot(ax, theta_old(idx_gt1), theta_new(idx_gt1)-1, ...
        LineStyle=plot_ls, Color=plot_colour, LineWidth=4.0, ...
        HandleVisibility='off');
end

% Plot diagonal line
plot(ax, [0, 1], [0, 1], LineStyle='--', Color=[0, 0, 0, 0.75], ...
      HandleVisibility='off');

% Legend
leg = legend(ax, Location='south east', Interpreter='latex', FontSize=25);
leg.ItemTokenSize = [80, 1];

hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = 0.0 : 0.25 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.125 : 0.25 : 1.0;

ax.YAxis.TickValues = 0.0 : 0.25 : 1.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.125 : 0.25 : 1.0;

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, 1.0];

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\theta_{\mathrm{old}}$';
ax.YAxis.Label.String = '$\theta_{\mathrm{new}}$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');

% Grid lines
% ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

% Filename
exportgraphics(fig, './PTC.pdf', ContentType='vector');

%-------------------------------------------------------------------------%
%%                               Functions                               %%
%-------------------------------------------------------------------------%
function read_write_data(filename_mat)
  % read_save_data(filename_mat)
  %
  % Reads the data from the solutions and saves to MATLAB data structure
  % thingy.

  if ~isfile(filename_mat)
    % .mat file does not exist so read data and save

    %-------------------%
    %     Read Data     %
    %-------------------%
    % Add main folder to path
    cd('../');

    % Folder name
    dir_data = sprintf('./data/%s/', 'run10_phase_reset_PTC_multi');
    % List all directories
    dirs = dir(dir_data);
    % Remove ./ and ../
    dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
    % Sub folder names
    dir_sub = {dirs.name};
    
    % Empty arrays for theta_old and theta_new
    theta_old_data = {};
    theta_new_data = {};
    % Empty array for A_perturb values
    A_perturb = [];
    theta_perturb = [];
    phi_perturb = [];
    
    % dir_sub_plot = dir_sub;
    dir_sub_plot = {dir_sub{1}; dir_sub{2}; dir_sub{9}; dir_sub{10}};
    
    % Cycle through each data directory
    for i = 1 : length(dir_sub_plot)
      % Sub folder name
      dir_read = {'run10_phase_reset_PTC_multi', dir_sub_plot{i}};
    
      % Bifurcation data
      bd_read = coco_bd_read(dir_read);
    
      % Read theta_old and theta_new
      theta_old_read = coco_bd_col(bd_read, 'theta_old');
      theta_new_read = coco_bd_col(bd_read, 'theta_new');
    
      % Read A_perturb
      A_perturb_read     = coco_bd_val(bd_read, 1, 'A_perturb');
    
      % Append to arrays
      theta_old_data{i} = theta_old_read;
      theta_new_data{i} = theta_new_read;
      A_perturb         = [A_perturb, A_perturb_read];
    
    end
    
    % Read directional vector components
    theta_perturb = coco_bd_val(bd_read, 1, 'theta_perturb');
    phi_perturb   = coco_bd_val(bd_read, 1, 'phi_perturb');
    % Directional vector
    d_vec = [cos(theta_perturb) * sin(phi_perturb);
             sin(theta_perturb) * sin(phi_perturb);
             cos(phi_perturb)];

    %--------------------%
    %     Write Data     %
    %--------------------%
    % Change back to abstract folder
    cd('./Nonlinear_Photonics_abstract/');

    % Save as matrix
    save(filename_mat, 'dir_sub_plot', 'theta_old_data', 'theta_new_data', ...
         'A_perturb', 'theta_perturb', 'phi_perturb', 'd_vec');
  
  end

end
