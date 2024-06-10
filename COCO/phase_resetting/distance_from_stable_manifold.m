%-------------------------------------------------------------------------%
%%        Calculate Distance and Perturbation from Stable Manifold       %%
%-------------------------------------------------------------------------%
% Load data
load('./data/initial_PO.mat');

% Find all points along stable manifold and periodic orbit where the
% Q value is "roughly" the same.

% Q values of periodic orbit
Q_PO = xbp_PO(:, 2);
% Q values of stable manfiold
Q_W = W_q_stable(:, 2);

% Empty array for differences
Q_diff = zeros(length(Q_PO), length(Q_W));

% Cycle through periodic orbit and stable manifold to find minimum
for i = 1 : length(Q_PO)
  for j = 1 : length(Q_W)
    % Calculate difference
    Q_diff(i, j) = abs(Q_PO(i) - Q_W(j));

  end
end

% Calculate min
min_val = zeros(length(Q_PO), 1);
min_idx = zeros(length(Q_PO), 1);

% Cycle through differences
for i = 1 : length(Q_PO)
  % Calculate minimum
  [min_temp, idx_temp] = min(Q_diff(i, :));

  % Update arrays
  min_val(i) = min_temp;
  min_idx(i) = idx_temp;
end

% Empty arrays for perturbation size and angle
A_perturb = zeros(length(Q_PO), 1);
theta_perturb = zeros(length(Q_PO), 1);

% Cycle through all Peroidic orbit values and calculate norm of distance
% vector
for i = 1 : length(Q_PO)
  % Get periodic orbit point
  vec_PO = [xbp_PO(i, 1), xbp_PO(i, 3)];

  % Get stable manifold point
  vec_W  = [W_q_stable(min_idx(i), 1), W_q_stable(min_idx(i), 3)];

  % Calculate difference
  vec_diff = vec_PO - vec_W;
  vec_diff = -vec_diff;
  A_diff = norm(vec_diff);

  % Angle of displacement vector
  % theta_diff = atan2(vec_diff(2), vec_diff(1));
  theta_diff = mod(atan2(vec_diff(2), vec_diff(1)), 2*pi);
  
  % if theta_diff < 0
  %   theta_diff = theta_diff + 2 * pi;
  % end

  % Update array
  A_perturb(i)     = A_diff;
  theta_perturb(i) = theta_diff;

end

% Get perturbation vector
d_vec = A_perturb .* [cos(theta_perturb), sin(theta_perturb)];
% d_vec = d_vec * 0.5;

% Plotting vector
G_plot = [xbp_PO(:, 1), xbp_PO(:, 1) + d_vec(:, 1)];
Q_plot = [xbp_PO(:, 2), xbp_PO(:, 2)];
I_plot = [xbp_PO(:, 3), xbp_PO(:, 3) + d_vec(:, 2)];

% List of all theta_old values
theta_old = linspace(0.0, 1.0, length(Q_PO));

% Sort
[~, sort_idx] = sort(theta_perturb);

%-------------------------------------------------------------------------%
%%                         Plot the Figures My G                         %%
%-------------------------------------------------------------------------%
plot_3D_PO_manifold_and_spokes(G_plot, Q_plot, I_plot, theta_perturb);

plot_A_and_theta_perturb_against_theta_old(A_perturb, theta_perturb, theta_old)

plot_2D_thingie(A_perturb, theta_perturb, theta_old);

%-------------------------------------------------------------------------%
%%                          PLOTTING FUNCTIONS                           %%
%-------------------------------------------------------------------------%
function plot_3D_PO_manifold_and_spokes(G_plot, Q_plot, I_plot, theta_perturb)
  %--------------------------%
  %     Load Matrix Data     %
  %--------------------------%
  load('./data/initial_PO.mat');

  %-----------------------------------------------------------------------%
  %%      Plot: Periodic Orbit and Stable Manifold with Displacements    %%
  %-----------------------------------------------------------------------%
  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
  % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(1); fig.Name = 'Initial Periodic Orbits'; clf;
  fig.Units = 'inches'; fig.Position = [1, 1, 16, 12]; fig.PaperSize = [16, 12];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Cycle through all points and plot perturbation vector towards stable
  % manifold.
  for i = 1 : length(Q_plot(:, 1))
    if theta_perturb(i) <= 0.5*pi
      % Plot
      plot3(ax, G_plot(i, :), Q_plot(i, :), I_plot(i, :), Color=colours(4, :), ...
            HandleVisibility='off');
    end

  end

  % Plot initial periodic orbit
  plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

  % Plot equilibrium points: x_{+}
  plot3(ax, x_pos(1), x_pos(2), x_pos(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
        MarkerEdgeColor='b', DisplayName='$q$');

  % Plot equilibrium points: x_{-}
  plot3(ax, x_neg(1), x_neg(2), x_neg(3), ...
        LineStyle='none', ...
        Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
        MarkerEdgeColor='r', DisplayName='$p$');

  % Plot equilibrium points: x_{0}
  plot3(ax, x_0(1), x_0(2), x_0(3), ...
        LineStyle='none', ...
        Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
        MarkerEdgeColor='r', DisplayName='$o$');

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
  % ax.XAxis.TickDirection = 'in';
  % ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
  % ax.XAxis.MinorTick = 'on';
  % ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;

  % Y-Axis
  % ax.YAxis.TickDirection = 'in';
  % ax.YAxis.TickValues = 0.0 : 2.0 : 8.0;
  % ax.YAxis.MinorTick = 'on';
  % ax.YAxis.MinorTickValues = 1.0 : 2.0 : 8.0;

  % Z-Axis
  % ax.ZAxis.TickDirection = 'in';
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

  %---------------------%
  %     Save Figure     %
  %---------------------%
  exportgraphics(fig, './images/Initial_PO_with_d_peturb.pdf', ContentType='vector');

end

function plot_A_and_theta_perturb_against_theta_old(A_perturb, theta_perturb, theta_old)
  %-----------------------------------------------------------------------%
  %%      Plot: Periodic Orbit and Stable Manifold with Displacements    %%
  %-----------------------------------------------------------------------%
  % Split data into theta_perturb <= pi/2 and everything else
  plot_idx = theta_perturb <= (0.5 * pi);

  theta_perturb_sub = theta_perturb(plot_idx);
  theta_perturb_sup = theta_perturb(~plot_idx);

  A_perturb_sub = A_perturb(plot_idx);
  A_perturb_sup = A_perturb(~plot_idx);

  theta_old_sub = theta_old(plot_idx);
  theta_old_sup = theta_old(~plot_idx);

  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
  % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(2); fig.Name = 'Perturbation and theta_old'; clf;
  fig.Units = 'inches'; fig.Position = [2, 2, 16, 12]; fig.PaperSize = [16, 12];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plot: theta_perturb <= 0.5 pi
  plot3(ax, theta_perturb_sub, theta_old_sub, A_perturb_sub, ...
        Color=colours(1, :), LineStyle='-', ...
        DisplayName='$\theta_{\mathrm{perturb}} \leq \pi / 2$');

  % Plot: theta_perturb <= 0.5 pi
  plot3(ax, theta_perturb_sup, theta_old_sup, A_perturb_sup, ...
        Color=[colours(1, :), 0.5], LineStyle='--', ...
        DisplayName='$\theta_{\mathrm{perturb}} > \pi / 2$');

  % Legend
  legend(ax, Interpreter='latex');

  % Turn of axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickDirection = 'in';
  ax.XAxis.TickValues = 0.0 : 0.25 * pi : 2 * pi;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 0.125*pi : 0.25 * pi : 2 * pi;
  ax.XAxis.TickLabels = {'0', '$\pi / 4$', '$\pi / 2$', '$3 \pi / 4$', ...
                         '$\pi$', '$5 \pi / 4$', '$3 \pi / 2$', ...
                         '$7 \pi / 4$', '$2 \pi$'};

  % Y-Axis
  ax.YAxis.TickDirection = 'in';
  ax.YAxis.TickValues = 0.0 : 0.1 : 1.0;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = 0.05 : 0.1 : 18.0;

  % Z-Axis
  ax.ZAxis.TickDirection = 'in';
  ax.ZAxis.TickValues = 0.0 : 3.0 : 15.0;
  ax.ZAxis.MinorTick = 'on';
  ax.ZAxis.MinorTickValues = 1.5 : 3.0 : 15.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 2*pi];
  ax.YAxis.Limits = [0.0, 1.0];
  ax.ZAxis.Limits = [0.0, 15.0];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$\theta_{\mathrm{perturb}}$';
  ax.ZAxis.Label.String = '$A_{\mathrm{perturb}}$';
  ax.YAxis.Label.String = '$\theta_{\mathrm{old}}$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  % ax.Title.String = 'Initial Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % 3D plot view
  view(45, 15.0);

  %---------------------%
  %     Save Figure     %
  %---------------------%
  exportgraphics(fig, './images/A_perturb_against_theta_old.pdf', ContentType='vector');

end

function plot_2D_thingie(A_perturb, theta_perturb, theta_old)
  %-----------------------------------------------------------------------%
  %%      Plot: Periodic Orbit and Stable Manifold with Displacements    %%
  %-----------------------------------------------------------------------%
  % Split data into theta_perturb <= pi/2 and everything else
  plot_idx = theta_perturb <= (0.5 * pi);

  theta_perturb_sub = theta_perturb(plot_idx);
  theta_perturb_sup = theta_perturb(~plot_idx);

  A_perturb_sub = A_perturb(plot_idx);
  A_perturb_sup = A_perturb(~plot_idx);

  theta_old_sub = theta_old(plot_idx);
  theta_old_sup = theta_old(~plot_idx);

  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
  % Default colour order (matplotlib)
  colours = colororder();

  fig = figure(3); fig.Name = 'Perturbation and theta_old'; clf;
  fig.Units = 'inches'; fig.Position = [3, 3, 10, 10]; fig.PaperSize = [10, 10];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Left axis
  yyaxis(ax, 'left');

  % Plot: theta_perturb <= 0.5 pi
  plot(ax, theta_old_sub, theta_perturb_sub, ...
       Color=colours(1, :), LineStyle='-', ...
       DisplayName='$\theta_{\mathrm{perturb}} \leq \pi / 2$');
  % Plot: theta_perturb <= 0.5 pi
  plot(ax, theta_old_sup, theta_perturb_sup, ...
       Color=[colours(1, :), 0.5], LineStyle='--', ...
       DisplayName='$\theta_{\mathrm{perturb}} > \pi / 2$');

  % Right axis
  yyaxis(ax, 'right');
  
  % Plot: theta_perturb <= 0.5 pi
  plot(ax, theta_old_sub, A_perturb_sub, ...
       Color=colours(2, :), LineStyle='-', ...
       HandleVisibility='off');
  % Plot: theta_perturb <= 0.5 pi
  plot(ax, theta_old_sup, A_perturb_sup, ...
       Color=[colours(2, :), 0.5], LineStyle='--', ...
       HandleVisibility='off');

  % Legend
  legend(ax, Interpreter='latex');

  % Turn of axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickDirection = 'in';
  ax.XAxis.TickValues = 0.0 : 0.1 : 1.0;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 0.05 : 0.1 : 18.0;

  % Y-Axis: theta_perturb
  ax.YAxis(1).TickDirection = 'in';
  ax.YAxis(1).TickValues = 0.0 : 0.25 * pi : 2 * pi;
  ax.YAxis(1).MinorTick = 'on';
  ax.YAxis(1).MinorTickValues = 0.125*pi : 0.25 * pi : 2 * pi;
  ax.YAxis(1).TickLabels = {'0', '$\pi / 4$', '$\pi / 2$', '$3 \pi / 4$', ...
                            '$\pi$', '$5 \pi / 4$', '$3 \pi / 2$', ...
                            '$7 \pi / 4$', '$2 \pi$'};

  % Y-Axis: A_perturb
  ax.YAxis(2).TickDirection = 'in';
  ax.YAxis(2).TickValues = 0.0 : 3.0 : 15.0;
  ax.YAxis(2).MinorTick = 'on';
  ax.YAxis(2).MinorTickValues = 1.5 : 3.0 : 15.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 1.0];

  % Y-Axis: theta_perturb
  ax.YAxis(1).Limits = [0.0, 2*pi];

  % Y-Axis: A_perturb
  ax.YAxis(2).Limits = [0.0, 15.0];

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$\theta_{\mathrm{old}}$';

  % Y-Axis: theta_perturb
  ax.YAxis(1).Label.String = '$\theta_{\mathrm{perturb}}$';

  % Y-Axis: A_perturb
  ax.YAxis(2).Label.String = '$A_{\mathrm{perturb}}$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  % ax.Title.String = 'Initial Periodic Orbit';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  exportgraphics(fig, './images/A_perturb_against_theta_old_2D.pdf', ContentType='vector');

end