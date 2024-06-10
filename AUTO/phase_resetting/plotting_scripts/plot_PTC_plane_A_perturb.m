%-------------------%
%     Read Data     %
%-------------------%
% Load data from .mat
load('../data/PTC_scan_data.mat')

theta_old_data = [];
theta_new_data = [];
A_perturb_data = [];

% Cycle through each data directory
for i = 1 : length(A_perturb)
  fprintf('%d \n', i);

  % Read theta_old and theta_new
  theta_old_read = theta_old{i};
  theta_old_read = theta_old_read';
  % theta_old_read = [theta_old_read(theta_old_read >= 1.0) - 1.0;
  %                   theta_old_read(theta_old_read < 1.0)];

  theta_new_read = theta_new{i};
  theta_new_read = theta_new_read';
  theta_new_read = [theta_new_read(theta_old_read >= 1.0);
                    theta_new_read(theta_old_read < 1.0)];

  % Read A_perturb
  A_perturb_read = A_perturb(i);
  A_perturb_read2 = A_perturb_read * ones(1, length(theta_new_read));

  % Pad arrays with NaN's so they're always the same length
  len_data = 2000 - length(theta_old_read);
  theta_old_read = padarray(theta_old_read, [0, len_data], NaN, 'pre');
  len_data = 2000 - length(theta_new_read);
  theta_new_read = padarray(theta_new_read, [0, len_data], NaN, 'pre');
  len_data = 2000 - length(A_perturb_read2);
  A_perturb_read2 = padarray(A_perturb_read2, [0, len_data], NaN, 'pre');

  % Append to arrays
  theta_old_data = [theta_old_data, theta_old_read'];
  theta_new_data = [theta_new_data, theta_new_read'];
  A_perturb_data = [A_perturb_data, A_perturb_read2'];

end

% size(theta_old_data)
% size(theta_new_data)
% size(A_perturb_data)

% % Read directional vector components
% theta_perturb = coco_bd_val(bd_read, 1, 'theta_perturb');
% phi_perturb   = coco_bd_val(bd_read, 1, 'phi_perturb');
% % Directional vector
% d_vec = [cos(theta_perturb) * sin(phi_perturb);
%           sin(theta_perturb) * sin(phi_perturb);
%           cos(phi_perturb)];

% %                             Plot Data                                 %
%-----------------------------------------------------------------------%
% Plotting colours
colours = colororder();

fig = figure(5); clf;
fig.Name = 'PTC Scans';
fig.Units = 'inches'; fig.Position = [0, 0, 16, 8]; fig.PaperSize = [16, 8];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 11;

%--------------%
%     Plot     %
%--------------%
hold(ax, 'on');

% % Cycle through and plot3
% for i = 1 : length(dir_sub)
%   % Grab data
%   theta_old_plot = theta_old_data{i};
%   theta_new_plot = theta_new_data{i};
%   A_perturb_plot = A_perturb_data{i};
% 
%   % size(theta_old_plot)
%   % size(theta_new_plot)
%   % size(A_perturb_plot)
% 
%   % Plot
%   plot3(ax, theta_old_plot, A_perturb_plot, theta_new_plot, ...
%         Color=colours(1, :), LineStyle='-');
% end

% Plot Surf plot (woah!)
surf(ax, theta_old_data, A_perturb_data, theta_new_data, ...
      FaceColor=colours(1, :), FaceAlpha=1, ...
      MeshStyle='column', LineStyle='-', EdgeColor=colours(2, :), ...
      LineWidth=0.5);

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%


%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\theta_{\mathrm{old}}$';
ax.YAxis.Label.String = '$A_{\mathrm{perturb}}$';
ax.ZAxis.Label.String = '$\theta_{\mathrm{new}}$';

%--------------------%
%     Axis Title     %
%--------------------%
% title_str = sprintf('Phase Transition Curve (PTC) with $\\vec{d} = (%.0f, %.0f, %.0f)$', d_vec(1), d_vec(2), d_vec(3));
% ax.Title.String = title_str;

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');
view(45, 15);

%---------------------%
%     Save Figure     %
%---------------------%
% if save_figure == true
%   % Filename
%   exportgraphics(fig, './images/PTC_plane_A_perturb.pdf', ContentType='vector');
% end