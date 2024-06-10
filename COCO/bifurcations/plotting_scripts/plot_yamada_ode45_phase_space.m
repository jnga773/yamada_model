% Clear plots
close all;

% Clear workspace
% clear;

% Add equation/functions to path
addpath('./functions/');

%-------------------------------------------------------------------------%
%%                     PLOT PHASE SPACE TRAJECTORIES                     %%
%-------------------------------------------------------------------------%
%---------------------------%
%     Plot: Phase Space     %
%---------------------------%
fig = figure(Name="ode45 Solution (Phase)", Units='inches', ...
             Position=[0, 0, 8, 8], PaperSize=[8, 8]);
ax = gca();

% Plot solutions
hold(ax, 'on');
for j = 1:length(x0_values)
    for k = 1:length(x0_values)
        for l = 1:length(x0_values)
            % Reshape data
            G = reshape(G_data(j, k, l, :), [1, length(t)]);
            Q = reshape(Q_data(j, k, l, :), [1, length(t)]);
            I = reshape(I_data(j, k, l, :), [1, length(t)]);
    
            % Plot
            plot(ax, G(:), I(:), LineStyle='-', LineWidth=0.5, ...
                 Color=[0.1216, 0.4667, 0.7059]);
            % plot(ax, t, u1(:), LineStyle='-', LineWidth=0.5, ...
            %      Color=[0.1216, 0.4667, 0.7059])
            % plot3(ax, G(:), Q(:), I(:), Color=[0.1216, 0.4667, 0.7059]);
        end
    end
end
hold(ax, 'off');

% Axis labels
ax.XAxis.Label.String = ['$G(t)$'];
ax.YAxis.Label.String = ['$I(t)$'];
% ax.ZAxis.Label.String = ['$I(t)$'];

% Axis title
ax.Title.String = ['Yamada Model'];

% Axis limits
% ax.XAxis.Limits = [-0.05, 1.05];
% ax.YAxis.Limits = [-0.05, 1.35];

% Tick params
ax.XAxis.TickDirection = 'in'; ax.YAxis.TickDirection = 'in';

% Figure stuff
grid(ax, 'on'); box(ax, 'on');
ax.GridLineWidth = 0.5;
ax.GridColor = 'black';
ax.GridAlpha = 0.25;

% Axis title
title_str = sprintf('Yamada Model with $\\left( \\gamma = %.2f, A = %.2f, B = %.2f, a = %.2f \\right)$', gamma, A, B, a);
ax.Title.String = title_str;

% view(45, 15.0)

% sgtitle('ode45 (Phase)');

% % Figure stuff
% exportgraphics(fig, './images/ode45_phase_profile.png', Resolution=600);