%-------------------------------------------------------------------------%
%%                               READ DATA                               %%
%-------------------------------------------------------------------------%
% Load matrix data
load('../data/level_set_scan.mat');

% Number of scans
num_scans = length(A_perturb);



%-------------------------------------------------------------------------%
%%                             PLOT SURFACE                              %%
%-------------------------------------------------------------------------%
% Plotting colours
colours = colororder();

% Create figure
fig = figure(5); clf;
fig.Name = 'Level Set Scan';
fig.Units = 'inches'; fig.Position = [0, 0, 16, 8]; fig.PaperSize = [16, 8];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 11;

%--------------%
%     Plot     %
%--------------%
hold(ax, 'on');

% Cycle through and plot3
for i = 1 : num_scans
  % Grab data
  theta_perturb_plot = theta_perturb{i};
  A_perturb_plot     = A_perturb{i};
  theta_new_plot     = theta_new(i) * ones(length(A_perturb_plot));

  % size(theta_old_plot)
  % size(theta_new_plot)
  % size(A_perturb_plot)

  % Plot
  plot3(ax, theta_perturb_plot, A_perturb_plot, theta_new_plot, ...
        Color=colours(1, :), LineStyle='-');
end

% % Plot Surf plot (woah!)
% surf(ax, theta_old_data, A_perturb_data, theta_new_data, ...
%       FaceColor=colours(1, :), FaceAlpha=1, ...
%       MeshStyle='column', LineStyle='-', EdgeColor=colours(2, :), ...
%       LineWidth=0.5);

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%


%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\theta_{\mathrm{perturb}}$';
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