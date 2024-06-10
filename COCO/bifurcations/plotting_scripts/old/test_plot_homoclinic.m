% Clear plots
% close('all');

%-----------------------------------------%
%     Plot High Period Periodic Orbit     %
%-----------------------------------------%
fig1 = figure(Name='High Period Periodic Orbit', Units='inches', ...
              Position=[0, 0, 12, 8], PaperSize=[12, 8]);
ax = gca();
ax.FontSize = 14.0; 

% Turn on axis hold
hold(ax, 'on');

plot(sol.tbp, sol.xbp(:, 1), LineStyle='-')
plot(sol.tbp, sol.xbp(:, 2), LineStyle='--')
plot(sol.tbp, sol.xbp(:, 3), LineStyle='-.')

% Turn off axis hold
hold(ax, 'off');

% Axis Labels
ax.XAxis.Label.String = ['$t$'];
ax.YAxis.Label.String = ['$x_{i}$'];

% Axis title
ax.Title.String = ['High Period Periodic Orbit'];

% Tick params
ax.XAxis.TickDirection = 'in'; ax.YAxis.TickDirection = 'in';

% Figure stuff
box(ax, 'on');
grid(ax, 'on');
ax.GridLineWidth = 0.5;
ax.GridColor = 'black';
ax.GridAlpha = 0.25;

%-------------------------------------%
%     Plot Sample Periodic Orbits     %
%-------------------------------------%
fig2 = figure(Name='Sample Periodic Orbits', Units='inches', ...
              Position=[0, 0, 12, 8], PaperSize=[12, 8]);
ax = gca();
ax.FontSize = 14.0;

hold(ax, 'on');
coco_plot_sol(run_old, [1:1:10, 15:5:100], '', 'x', 1, 'x', 2)
hold(ax, 'off');

% Axis Labels
ax.XAxis.Label.String = ['$G$'];
ax.YAxis.Label.String = ['$Q$'];

% Axis title
ax.Title.String = ['High Period Periodic Orbit'];

% Tick params
ax.XAxis.TickDirection = 'in'; ax.YAxis.TickDirection = 'in';

% Figure stuff
box(ax, 'on');
grid(ax, 'on');
ax.GridLineWidth = 0.5;
ax.GridColor = 'black';
ax.GridAlpha = 0.25;