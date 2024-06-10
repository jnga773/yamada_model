% Close all plots
% close('all');

%--------------------------------------%
%     Plot: Panel (f) New Solution     %
%--------------------------------------%
fig = figure(6); clf;
fig.Name = '(f) New Solution';
ax = gca();

hold(ax, 'on');
% Plot COCO solution: 8
coco_plot_sol(run2, 8, '', 'x', 2, 'x', 3)
coco_plot_sol(run3, 1, '', 'x', 2, 'x', 3)
hold(ax, 'off');

% Axis Labels
ax.XAxis.Label.String = ['$x_{2}$'];
ax.YAxis.Label.String = ['$x_{3}$'];

% Axis title
ax.Title.String = ['(f) New Solution'];

% Tick params
ax.XAxis.TickDirection = 'in'; ax.YAxis.TickDirection = 'in';

% Figure stuff
box(ax, 'on');
grid(ax, 'on');
ax.GridLineWidth = 0.5;
ax.GridColor = 'black';
ax.GridAlpha = 0.25;

exportgraphics(fig, './images/(f) New Solution.png', Resolution=800);
