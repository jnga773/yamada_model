% Filename of .nat file
% mat_file = '../data_mat/PTC_scan.mat';
mat_file = '../data_mat/PTC_scan_G.mat';
% mat_file = '../data_mat/PTC_scan_I.mat';

% Filename of figure
filename_out = '../images/PTC_scan.pdf';
% filename_out = '../images/PTC_scan_I_small_range.pdf';

% Save figure or nah?
save_figure = false;
% save_figure = true;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
% Load data from .mat
load(mat_file, 'theta_perturb');

% Read directional vector components
% Directional vector
d_vec = [cos(theta_perturb); 0; sin(theta_perturb)];

%------------------%
%     Pad Data     %
%------------------%
% [theta_old_data, theta_new_data, A_perturb_data] = read_data(mat_file);
load(mat_file);

% [theta_old_data, theta_new_data, A_perturb_data] = pad_data(theta_old_data, theta_new_data, A_perturb_data);

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
% Plot theta_old > 1
for i = 1 : length(A_perturb_gt1(1, :))
  % Grab data
  theta_old_plot = theta_old_gt1(:, i);
  theta_new_plot = theta_new_gt1(:, i);
  A_perturb_plot = A_perturb_gt1(:, i);

  % Plot
  plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
        Color=colours(1, :), LineStyle='-');
end

% Plot theta_old < 1
for i = 1 : length(A_perturb_lt1(1, :))
  % Grab data
  theta_old_plot = theta_old_lt1(:, i);
  theta_new_plot = theta_new_lt1(:, i);
  A_perturb_plot = A_perturb_lt1(:, i);

  % Plot
  plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
        Color=colours(2, :), LineStyle='-');
end

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% % Plot Surf plot (woah!)
% [X_gt1, Y_gt1, Z_gt1] = matrix_data(theta_old_gt1, theta_new_gt1, A_perturb_gt1);
% 
% surf(ax, X_gt1, Y_gt1, Z_gt1, ...
%       FaceColor=colours(1, :), FaceAlpha=1, ...
%       MeshStyle='column');
% 
% % Plot Surf plot (woah!)
% [X_lt1, Y_lt1, Z_lt1] = matrix_data(theta_old_lt1, theta_new_lt1, A_perturb_lt1);
% 
% surf(ax, X_lt1, Y_lt1, Z_lt1, ...
%       FaceColor=colours(1, :), FaceAlpha=1, ...
%       MeshStyle='column');
% 
% % Plot Surf plot (woah!)
% [X_gt1, Y_gt1, Z_gt1] = matrix_data(theta_old_gt1, theta_new_gt1-1.0, A_perturb_gt1);
% 
% surf(ax, X_gt1, Y_gt1, Z_gt1, ...
%       FaceColor=colours(1, :), FaceAlpha=1, ...
%       MeshStyle='column');
% 
% % Plot Surf plot (woah!)
% [X_lt1, Y_lt1, Z_lt1] = matrix_data(theta_old_lt1, theta_new_lt1-1.0, A_perturb_lt1);
% 
% surf(ax, X_lt1, Y_lt1, Z_lt1, ...
%       FaceColor=colours(1, :), FaceAlpha=1, ...
%       MeshStyle='column');
% 
% % Shading of surface
% shading(ax, 'interp');

% Colorbar
cbar = colorbar();

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
% ax.YAxis.Limits = [0.0, 2.0];
% ax.ZAxis.Limits = [0.0, 3.0];

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\theta_{\mathrm{old}}$';
ax.YAxis.Label.String = '$A_{\mathrm{perturb}}$';
ax.ZAxis.Label.String = '$\theta_{\mathrm{new}}$';

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
view(-45, 15);

%---------------------%
%     Save Figure     %
%---------------------%
if save_figure == true
  % Filename
  exportgraphics(fig, filename_out, ContentType='vector');
end

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
function [theta_old_out, theta_new_out, A_perturb_out] = pad_data(theta_old_in, theta_new_in, A_perturb_in)
  % theta_old_out, theta_new_out, A_perturb_out = pad_data(theta_old_in, theta_new_in, A_perturb_in)
  %
  % Pads the data arrays with NaNs

  % Find lengths of all data arrays
  array_lengths = [];
  for i = 1 : length(theta_old_in)
    array_lengths = [array_lengths; length(theta_old_in{i})];
  end

  % Max data length
  max_length = max(array_lengths);

  % Empty arrays for output data
  theta_old_out = {};
  theta_new_out = {};
  A_perturb_out = {};

  % Cycle through and pad data
  for i = 1 : length(theta_old_in)
    % Length to pad
    pad_length = max_length - length(theta_old_in{i});

    A_read = A_perturb_in{i};
    old_read = theta_old_in{i}';
    new_read = theta_new_in{i}';

    % Pad arrays
    theta_old_out{i} = padarray(old_read, [pad_length, 0], old_read(1), 'pre');
    theta_new_out{i} = padarray(new_read, [pad_length, 0], new_read(1), 'pre');
    A_perturb_out{i} = padarray(A_read  , [pad_length, 0], A_read(1)  , 'pre');

  end

  theta_old_out = theta_old_out';
  theta_new_out = theta_new_out';
  A_perturb_out = A_perturb_out';

end

function [theta_old_out, theta_new_out, A_perturb_out] = read_data(mat_file_in)
  % [theta_old_out, theta_new_out, A_perturb_out] = read_data()
  %
  % Reads and tidies up the data from the Python generated .mat file.

  load(mat_file_in, 'theta_old_gt1', 'theta_old_lt1', ...
                    'theta_new_gt1', 'theta_new_lt1', ...
                    'A_perturb');

  % Empty arrays for output data
  theta_old_out = {};
  theta_new_out = {};
  A_perturb_out = {};

  % Cycle through each data directory
  for i = 1 : length(A_perturb)

    % Read theta_old and theta_new
    theta_old_lt1_read = theta_old_lt1{i};
    theta_old_gt1_read = theta_old_gt1{i};

    theta_new_lt1_read = theta_new_lt1{i};
    theta_new_gt1_read = theta_new_gt1{i};

    % Check end and start values of theta_new_array
    lt1 = [theta_new_lt1_read(1), theta_new_lt1_read(end)];
    gt1 = [theta_new_gt1_read(1), theta_new_gt1_read(end)];

    mean_diff_1 = mean(abs(lt1(1) - gt1(1)));
    mean_diff_end = mean(abs(lt1(end) - gt1(end)));
    % Mean these two mfers
    mean_diff = 0.5 * (mean_diff_1 + mean_diff_end);

    if mean_diff < 1e-3
      % Difference is pretty small
      % Interpolate / mix data up
      theta_old_read = [theta_old_gt1_read, theta_old_lt1_read];
      theta_new_read = [theta_new_gt1_read, theta_new_lt1_read];

      % Sort by theta_old
      [theta_old_read, sort_idx] = sort(theta_old_read);
      theta_new_read = theta_new_read(sort_idx);

      % Read A_perturb
      A_perturb_read = A_perturb(i) * ones(length(theta_old_read), 1);
    else
      % If there is a gap between the gt1 and lt1 data, append them with a
      % NaN in the centre.
      theta_old_read = [theta_old_gt1_read, NaN, theta_old_lt1_read]';
      theta_new_read = [theta_new_gt1_read, NaN, theta_new_lt1_read]';

      % Read A_perturb
      A_perturb_read = A_perturb(i) * ones(length(theta_old_read), 1);
    end

    % Append to arrays 
    theta_old_out{i} = theta_old_read;
    theta_new_out{i} = theta_new_read;
    A_perturb_out{i} = A_perturb_read;
  end

  % Pad data
  [theta_old_out, theta_new_out, A_perturb_out] = pad_data(theta_old_out, theta_new_out, A_perturb_out);

end

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