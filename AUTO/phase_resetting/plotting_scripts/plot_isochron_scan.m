close all; clear all; clc;

%---------------------------------------%
%     Save Data to MATLAB .mat file     %
%---------------------------------------%
% Save data
load('../data_mat/isochron_scan.mat', ...
    'iso1_data', 'iso2_data', 'iso3_data', ...
    'parameters');

%--------------------------------------%
%     Read Data: Unperturbed Orbit     %
%--------------------------------------%
% Read unperturbed periodic orbit data
load('./initial_PO_MATLAB.mat');

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

% % Plot unperturbed periodic orbit
% plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
%     LineStyle='-', Color=colours(3, :), ...
%     DisplayName='$\Gamma$');
% 
% % Plot stable manifold of q / x_{+}
% plot3(ax, W_q_stable(:, 1), W_q_stable(:, 2), W_q_stable(:, 3), ...
%     Color=colours(1, :), ...
%     DisplayName='$W^{s}(p)$');

% Cycle through data and plot
for i = 1 : length(iso1_data)

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
ax.XAxis.Limits = [-2.0, 6.0];
ax.YAxis.Limits = [-2.0, 6.0];
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
