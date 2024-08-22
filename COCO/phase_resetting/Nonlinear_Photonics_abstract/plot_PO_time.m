clear all; close all; clc;

%---------------------%
%%     Read Data     %%
%---------------------%
% Check to see if 'data_PO_time.mat' exists. If not, read data and create
% it
filename_data = './data_PO_time.mat';

read_write_data(filename_data);

% Load data
load(filename_data);

%-------------------------------------------------------------------------%
%%          Nonlinear Photonics 2024 Abstract - Periodic Orbit           %%
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Initial Periodic Orbits (Time)';
fig.Units = 'inches'; fig.Position = [3, 3, 16, 4]; fig.PaperSize = [16, 4];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 18;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plot original periodic orbit
plot(ax, t_PO_plot, x_PO_plot, Color=colours(3, :), ...
     LineWidth=2.0, DisplayName='$\Gamma$');

% Plot segment 4
plot(ax, t4, x4(:, 3), Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=2.0, DisplayName='Segment 4');

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 500.0];
ax.YAxis.Limits = [0.0, 21];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 100.0 : 500.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 50.0 : 100.0 : 500.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 5.0 : 20.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 2.5 : 5.0 : 20.0;

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$t$';
ax.YAxis.Label.String = '$I$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% Grid lines
ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

exportgraphics(fig, './periodic_orbit_time.pdf', ContentType='vector');

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

    % Base periodic orbit solution
    sol_PO = coll_read_solution('initial_PO', 'run06_initial_periodic_orbit', 1);
    x_PO   = sol_PO.xbp;
    t_PO   = sol_PO.tbp;
    
    % append x_PO and t_PO 25 times
    x_PO_plot = [];
    t_PO_plot = [];
    
    t_max = 0.0;
    
    for i = 1 : 25
      % Append x_PO
      x_PO_plot = [x_PO_plot; x_PO(:, 3)];
    
      % Append t_PO
      t_PO_plot = [t_PO_plot; t_max + t_PO];
    
      t_max = t_PO_plot(end);
    end
    
    % Phase reset orbit
    sol4 = coll_read_solution('seg4', 'run09_phase_reset_perturbation', 22);
    x4   = sol4.xbp;
    t4   = sol4.tbp;
    t4 = 25 * t4;

    %--------------------%
    %     Write Data     %
    %--------------------%
    % Change back to abstract folder
    cd('./Nonlinear_Photonics_abstract/');

    % Save as matrix
    save(filename_mat, 'x_PO_plot', 't_PO_plot', 'x4', 't4');
  
  end

end