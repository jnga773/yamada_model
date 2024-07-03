% Clear stuff
clear all;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
mat_file = '../data_mat/PTC_scan_G.mat';
% Load data from .mat
load(mat_file);

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
% Plotting colours
colours = colororder();

fig = figure(5); clf;
fig.Name = 'PTC Scans';
fig.Units = 'inches'; fig.Position = [0, 0, 16, 8]; fig.PaperSize = [16, 8];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 11;

hold(ax, 'on');

%-----------------------%
%     Plot: 3D Plot     %
%-----------------------%
% % Plot: Hole (theta_old < 1)
% data_plot = data_hole_lt1;
% for i = 1 : length(data_plot.A_perturb(1, :))
%   % Grab data
%   theta_old_plot = data_plot.theta_old(:, i);
%   theta_new_plot = data_plot.theta_new(:, i);
%   A_perturb_plot = data_plot.A_perturb(:, i);
% 
%   % Plot
%   plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
%         Color=colours(1, :), LineStyle='-');
% end
% 
% % Plot: Hole (theta_old > 1)
% data_plot = data_hole_gt1;
% for i = 1 : length(data_plot.A_perturb(1, :))
%   % Grab data
%   theta_old_plot = data_plot.theta_old(:, i);
%   theta_new_plot = data_plot.theta_new(:, i);
%   A_perturb_plot = data_plot.A_perturb(:, i);
% 
%   % Plot
%   plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
%         Color=colours(1, :), LineStyle='-');
% end
% 
% % Plot: Before hole
% data_plot = data_before_hole;
% for i = 1 : length(data_plot.A_perturb(1, :))
%   % Grab data
%   theta_old_plot = data_plot.theta_old(:, i);
%   theta_new_plot = data_plot.theta_new(:, i);
%   A_perturb_plot = data_plot.A_perturb(:, i);
% 
%   % Plot
%   plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
%         Color=colours(1, :), LineStyle='-');
% end
% 
% % Plot: Before hole
% data_plot = data_after_hole;
% for i = 1 : length(data_plot.A_perturb(1, :))
%   % Grab data
%   theta_old_plot = data_plot.theta_old(:, i);
%   theta_new_plot = data_plot.theta_new(:, i);
%   A_perturb_plot = data_plot.A_perturb(:, i);
% 
%   % Plot
%   plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
%         Color=colours(1, :), LineStyle='-');
% end

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = matrix_data(data_hole_lt1.theta_old, data_hole_lt1.theta_new, data_hole_lt1.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: Hole (theta_old > 1)
[X, Y, Z] = matrix_data(data_hole_gt1.theta_old, data_hole_gt1.theta_new, data_hole_gt1.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: Before hole
[X, Y, Z] = matrix_data(data_before_hole.theta_old, data_before_hole.theta_new, data_before_hole.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: After hole
[X, Y, Z] = matrix_data(data_after_hole.theta_old, data_after_hole.theta_new, data_after_hole.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

%-----------------------------------------%
%     Plot: Surface (One Level Lower)     %
%-----------------------------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = matrix_data(data_hole_lt1.theta_old, data_hole_lt1.theta_new - 1.0, data_hole_lt1.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: Hole (theta_old > 1)
[X, Y, Z] = matrix_data(data_hole_gt1.theta_old, data_hole_gt1.theta_new - 1.0, data_hole_gt1.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: Before hole
[X, Y, Z] = matrix_data(data_before_hole.theta_old, data_before_hole.theta_new - 1.0, data_before_hole.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: After hole
[X, Y, Z] = matrix_data(data_after_hole.theta_old, data_after_hole.theta_new - 1.0, data_after_hole.A_perturb);
surf(ax, X, Y, Z, MeshStyle='column');

%---------------------------%
%     Surface: Settings     %
%---------------------------%
% Shading of surface
shading(ax, 'interp');

% Colorbar
cbar = colorbar();

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, max(A_perturb)];
ax.ZAxis.Limits = [-1.0, 3.0];

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\theta_{\mathrm{old}}$';
ax.YAxis.Label.String = '$A_{\mathrm{perturb}}$';
ax.ZAxis.Label.String = '$\theta_{\mathrm{new}}$';

%--------------------%
%     Axis Title     %
%--------------------%
d_vec = [cos(theta_perturb); 0.0; sin(theta_perturb)];
title_str = sprintf('Phase Transition Curve (PTC) with $\\vec{d} = (%.0f, %.0f, %.0f)$', d_vec(1), d_vec(2), d_vec(3));
ax.Title.String = title_str;

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');
view(-45, 15);

%---------------------%
%     Save Figure     %
%---------------------%
% if save_figure == true
%   % Filename
%   exportgraphics(fig, filename_out, ContentType='vector');
% end

%-------------------------------------------------------------------------%
%%                           Data Functions                              %%
%-------------------------------------------------------------------------%
function [theta_old, A_perturb, theta_new] = matrix_data(theta_old_in, theta_new_in, A_perturb_in)
  % [theta_old, theta_new, A_perturb] = matrix_data(theta_old_in, theta_new_in, A_perturb_in)
  %
  % Turns the cells of arrays into a matrix

  theta_old = [];
  theta_new = [];
  A_perturb = [];

  % Cycle through cells and append as columns to matrix
  size_array = size(A_perturb_in);

  for i = 1 : size_array(1)
    theta_old_read = theta_old_in(i, :);
    theta_new_read = theta_new_in(i, :);
    A_perturb_read = A_perturb_in(i, :);

    theta_old = [theta_old, theta_old_read'];
    theta_new = [theta_new, theta_new_read'];
    A_perturb = [A_perturb, A_perturb_read'];

  end

end