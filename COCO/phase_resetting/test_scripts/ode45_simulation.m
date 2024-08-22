%=========================================================================%
%                     YAMADA MODEL (ode45 Simulation)                     %
%=========================================================================%
% We compute the phase resetting of an attracting periodic orbit of the
% Yamada model:
%                     G' = \gamma (A - G - G I) ,
%                     Q' = \gamma (B - Q - a Q I) ,
%                     I' = (G - Q - 1) I ,
% where G is the gain, Q is the absorption, and I is the intensity of the
% laser. The system is dependent on four parameters: the pump current on
% the gain, A (or A); the relative absoprtion, B and a; and the decay
% time of the gain, \gamma.

% Clear plots
close('all');

% Clear workspace
clear;
clc;

% Add functions to path
addpath('./functions/fields/hardcoded/');

%-------------------%
%     Read Data     %
%-------------------%
% Load in data from the ./data_mat/initial_PO.mat file
load('./data_mat/initial_PO.mat');

% Normalise orbit data
T_PO = max(tbp_PO);
tbp_norm = tbp_PO / T_PO;

% Pick theta_old point
theta_old = 0.6;
% Perturbation amplitude
A_perturb = 5.0;
% Pertubation direction
d_perturb = [0.0, 0.0, 1.0];
% Periodicity
k         = 500;

% Find point along periodic orbit closest to theta_old
[~, min_idx] = min(abs(tbp_norm - theta_old));

% Take this point as the initial condition
x_theta_old = xbp_PO(min_idx, :);

% Perturb to get initial condition
x_init = x_theta_old + (A_perturb * d_perturb);

%---------------------%
%     Solve ode45     %
%---------------------%
% Temporary function to solve
func_temp = @(t_in, x_in) yamada(x_in, p_PO);

% Time span
t_span = [0, k * T_PO];

% Solve that shiet
[t_solve, x_solve] = ode45(func_temp, t_span, x_init);
t_solve = t_solve / T_PO;
t_solve = t_solve + theta_old;

%-------------------------------------------------------------------------%
%%                                Plot: 3D                               %%
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

fig = figure(1); fig.Name = 'Initial Periodic Orbits'; clf;
fig.Units = 'inches'; fig.Position = [3, 3, 8, 8]; fig.PaperSize = [8, 8];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plotting colours
colours = colororder();

% Plot ode45 solution
plot3(ax, x_solve(:, 1), x_solve(:, 2), x_solve(:, 3), ...
      LineStyle='-', Color='k', ...
      DisplayName='ode45');

% Plot initial periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
      LineStyle='-', Color=colours(3, :), ...
      DisplayName='$\Gamma$');

% Plot equilibrium points: x_{+}
plot3(ax, x_pos(1), x_pos(2), x_pos(3), ...
      LineStyle='none', ...
      Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
      MarkerEdgeColor='b', DisplayName='$q$');

% % Plot equilibrium points: x_{-}
% plot3(ax, x_neg(1), x_neg(2), x_neg(3), ...
%       LineStyle='none', ...
%       Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
%       MarkerEdgeColor='r', DisplayName='$p$');

% % Plot equilibrium points: x_{0}
% plot3(ax, x_0(1), x_0(2), x_0(3), ...
%       LineStyle='none', ...
%       Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
%       MarkerEdgeColor='r', DisplayName='$o$');

% Plot stable manifold of q / x_{+}
plot3(ax, W_q_stable(:, 1), W_q_stable(:, 2), W_q_stable(:, 3), ...
      Color=colours(1, :), ...
      DisplayName='$W^{s}(p)$');

% Legend
legend(ax, Interpreter='latex');

% Turn of axis hold
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 1.0 : 5.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 1.0 : 5.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

% Z-Axis
ax.ZAxis.TickDirection = 'in';
% ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
% ax.ZAxis.MinorTick = 'on';
% ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 5.0];
ax.YAxis.Limits = [0.0, 4.0];
ax.ZAxis.Limits = [0.0, ceil(max(xbp_PO(:, 3)))];

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

%-------------------------------------------------------------------------%
%%                                Plot: 2D                               %%
%-------------------------------------------------------------------------%
%-----------------------%
%     Plotting Data     %
%-----------------------%
% Add k number of periods to data
xbp_PO_plot = [xbp_PO(1, :)];
tbp_PO_plot = [tbp_norm(1)];

t_max = 0.0;

for i = 1 : k
  % Append x_PO
  xbp_PO_plot = [xbp_PO_plot; xbp_PO(2:end, :)];
  % Append t_PO
  tbp_PO_plot = [tbp_PO_plot; t_max + tbp_norm(2:end)];
  t_max = tbp_PO_plot(end);
end

%----------------------%
%     Figure Setup     %
%----------------------%
% Default colour order (matplotlib)
colours = colororder();

fig = figure(3); fig.Name = 'Initial Periodic Orbits'; clf;
fig.Units = 'inches'; fig.Position = [3, 3, 16, 8]; fig.PaperSize = [16, 8];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plotting colours
colours = colororder();

% Plot ode45 solution
plot(ax, t_solve, x_solve(:, 3), ...
      LineStyle='-', Color=[0.0, 0.0, 0.0, 0.5], LineWidth=2.5, ...
      DisplayName='ode45');

% Plot initial periodic orbit
plot(ax, tbp_PO_plot, xbp_PO_plot(:, 3), ...
     LineStyle='-', Color=[colours(3, :), 0.35], LineWidth=1.5, ...
     DisplayName='$\Gamma$');

% Legend
legend(ax, Interpreter='latex');

% Turn of axis hold
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
% ax.XAxis.TickValues = 0.0 : 1.0 : 5.0;
% ax.XAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

% Y-Axis
% ax.YAxis.TickValues = 0.0 : 1.0 : 5.0;
% ax.YAxis.MinorTick = 'on';
% ax.YAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, max(t_solve)];
% ax.YAxis.Limits = [0.0, 4.0];

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$t / T_{\Gamma}$';
ax.YAxis.Label.String = '$I(t)$';

%--------------------%
%     Axis Title     %
%--------------------%
ax.Title.String = 'Initial Periodic Orbit';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');
