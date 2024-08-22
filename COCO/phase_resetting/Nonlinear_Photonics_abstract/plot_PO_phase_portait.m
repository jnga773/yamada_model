clear all; close all; clc;

%---------------------%
%%     Read Data     %%
%---------------------%
% Check to see if 'data_PO_phase.mat' exists. If not, read data and create
% it
filename_data = './data_PO_phase.mat';

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
fig.Name = 'Initial Periodic Orbits (Phase)';
fig.Units = 'inches'; fig.Position = [3, 3, 8, 7]; fig.PaperSize = [8, 7];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 18;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plot segment 4
plot3(ax, x4(:, 1), x4(:, 2), x4(:, 3), Color=[0.0, 0.0, 0.0, 0.5], ...
      LineWidth=1.0, DisplayName='Segment 4');

% Plot original periodic orbit
plot3(ax, x_PO(:, 1), x_PO(:, 2), x_PO(:, 3), Color=colours(3, :), ...
      LineWidth=4.0, DisplayName='$\Gamma$');

% Plot equilibrium point
plot3(ax, x_ep(1), x_ep(2), x_ep(3), Color='k', Marker='*', MarkerSize=20, ...
      MarkerFaceColor='k', MarkerEdgecolor='k', ...
      HandleVisibility='off');

% Plot stable manifold
plot3(ax, W_stable_p(:, 1), W_stable_p(:, 2), W_stable_p(:, 3), ...
      Color=colours(1, :), LineWidth=4.0, DisplayName='$W^{s}(p)$');

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.YAxis.Limits = [0.0, 4.0];
ax.ZAxis.Limits = [0.0, 21];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 2.0 : 6.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 1.0 : 2.0 : 6.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 1.0 : 4.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.5 : 1.0 : 4.0;

% Z-Axis
ax.ZAxis.TickDirection = 'in';
ax.ZAxis.TickValues = 0.0 : 5.0 : 20.0;
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = 2.5 : 5.0 : 20.0;

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$G$';
ax.YAxis.Label.String = '$Q$';
ax.ZAxis.Label.String = '$I$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% Grid lines
ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

% 3D plot view
view(45, 10.0);

exportgraphics(fig, './periodic_orbit_3D.pdf', ContentType='vector');

%-------------------------------------------------------------------------%
%%                               Functions                               %%
%-------------------------------------------------------------------------%
function x_out = calc_stable_manifold(x_in, p_in)
  % x_out = calc_stable_manifold()
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Add function folder to path
  addpath('./functions/');

  % Jacobian
  J_stable = yamada_DFDX(x_in, p_in);

  % Calculate eigenvalues and eigenvectors
  [eigvec, eigval] = eig(J_stable);

  % Stable eigenvector
  vec_s = eigvec(:, 3);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps1 = -0.01;
  % Time span
  t_span1 = [0.0, -16.5];

  % Initial vector
  x_init1 = x_in + (eps1 * vec_s);

  % Integrate using ode45
  [~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p_in), t_span1, x_init1);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps2 = 0.01;
  % Time span
  t_span2 = [0.0, -20.0];

  % Initial vector
  x_init2 = x_in + (eps2 * vec_s);

  % Integrate using ode45
  [~, W2] = ode45(@(t_in, x0_in) yamada(x0_in, p_in), t_span2, x_init2);

  %----------------%
  %     Output     %
  %----------------%
  x_out = [flip(W2); W1];

end

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
    
    % Phase reset orbit
    sol4 = coll_read_solution('seg4', 'run09_phase_reset_perturbation', 22);
    x4   = sol4.xbp;
    
    % Read equilibrium point data
    sol_EP = ep_read_solution('singularity', 'run09_phase_reset_perturbation', 22);
    
    % Equiliubrium point
    x_ep = sol_EP.x;
    
    %-----------------------------------%
    %     Calculate Stable Manifold     %
    %-----------------------------------%
    W_stable_p = calc_stable_manifold(x_ep, sol_EP.p);

    %--------------------%
    %     Write Data     %
    %--------------------%
    % Change back to abstract folder
    cd('./Nonlinear_Photonics_abstract/');

    % Save as matrix
    save(filename_mat, 'x_PO', 'x4', 'x_ep', 'W_stable_p');
  
  end

end