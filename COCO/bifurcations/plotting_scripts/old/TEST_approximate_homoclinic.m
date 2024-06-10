% Close all plots
% close('all');

%------------------------------------------------%
%     Plot: Panel (g) Approximate Homoclinic     %
%------------------------------------------------%
fig = figure(7); clf;
fig.Name = '(g) Approximate Homoclinic';
ax = gca();

hold(ax, 'on');
% Plot COCO solution: bifurcation diagrams
thm = struct('ustab', '', 'xlab', 'p_1', 'ylab', 'p_2');
thm.lspec = {'k', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 8};
coco_plot_bd(thm, run4, 'p1', 'p2')
hold(ax, 'off');

% Axis Labels
ax.XAxis.Label.String = ['$p_{1}$'];
ax.YAxis.Label.String = ['$p_{2}$'];

% Axis title
ax.Title.String = ['(g) Approximate Homoclinic'];

% Tick params
ax.XAxis.TickDirection = 'in'; ax.YAxis.TickDirection = 'in';

% Figure stuff
box(ax, 'on');
grid(ax, 'on');
ax.GridLineWidth = 0.5;
ax.GridColor = 'black';
ax.GridAlpha = 0.25;

exportgraphics(fig, './images/(g) Approximate Homoclinic.png', Resolution=800);
