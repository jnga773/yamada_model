%-------------------------------------------------------------------------%
%%                        Compute Stable Manifold                        %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Which run this continuation continues from
run_old = 'run06_initial_periodic_orbit';

% Continuation point
label_old = 1;

% Print to console
fprintf("~~~ Initial Periodic Orbit: (compute_stable_manifold.m) ~~~ \n");
fprintf('Compute stable manifold of equilibrium point \n');
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-----------------------%
%     Do Some Stuff     %
%-----------------------%
% plot_stable_manifold(run_old, label_old);

%-------------------%
%     Read Data     %
%-------------------%
% Read EP solution
sol_EP   = ep_read_solution('xpos', run_old, label_old);

% Equilibrium point
x_ss = sol_EP.x;
% Parameters
p    = sol_EP.p;

%------------------------------%
%     Calculate EigenStuff     %
%------------------------------%
% Jacobian
J_stable = yamada_DFDX(x_ss, p);

% Calculate eigenvalues and eigenvectors
[eigvec, eigval] = eig(J_stable);

% Indices for stable eigenvectors (eigval < 0)
stable_index = find(diag(eigval) < 0);

% Stable eigenvector
v1_s = eigvec(:, stable_index(1));
v2_s = eigvec(:, stable_index(2));

% Create a mesh of initial conditions
[s, t] = meshgrid(linspace(-0.1, 0.1, 20), linspace(-0.1, 0.1, 20));

%------------------------------------------------%
%     Calculate Things: Positive y direction     %
%------------------------------------------------%
% Small distance
eps1 = -0.01;
% Time span
t_span1 = [0.0, -16.0];

% Initial vector
% x_init1 = x_ss + (s(:) * v1_s) + (t(:) * v2_s);
x_init = x_ss + (s(:) * real(v1_s)) + (t(:) * imag(v1_s));

% Integrate using ode45
[~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span1, x_init1);

%-------------------------------------------------------------------------%
%%                           Some Functions Ow                           %%
%-------------------------------------------------------------------------%
function x_out = calc_stable_manifold(run_in, label_in)
  % x_out = calc_stable_manifold(run_in, label_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read EP solution
  sol_EP   = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium point
  x_ss = sol_EP.x;
  % Parameters
  p    = sol_EP.p;
  
  %------------------------------%
  %     Calculate EigenStuff     %
  %------------------------------%
  % Jacobian
  J_stable = yamada_DFDX(x_ss, p);

  % Calculate eigenvalues and eigenvectors
  [eigvec, eigval] = eig(J_stable);

  % Indices for stable eigenvectors (eigval < 0)
  stable_index = find(diag(eigval) < 0);

  % Stable eigenvector
  v1_s = eigvec(:, stable_index(1));
  v2_s = eigvec(:, stable_index(2));

  % Create a mesh of initial conditions
  [s, t] = meshgrid(linspace(-0.1, 0.1, 20), linspace(-0.1, 0.1, 20));

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps1 = -0.01;
  % Time span
  t_span1 = [0.0, -16.0];

  % Initial vector
  x_init1 = x_ss + (s(:) * v1_s) + (t(:) * v2_s);

  % Integrate using ode45
  [~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span1, x_init1);

  % %------------------------------------------------%
  % %     Calculate Things: Positive y direction     %
  % %------------------------------------------------%
  % % Small distance
  % eps2 = 0.01;
  % % Time span
  % t_span2 = [0.0, -20.0];
  % 
  % % Initial vector
  % x_init2 = x_ss + (eps2 * vec1_s);
  % 
  % % Integrate using ode45
  % [~, W2] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span2, x_init2);

  %----------------%
  %     Output     %
  %----------------%
  % x_out = [flip(W2); W1];
  x_out = W1;

end

% function write_data(W_stable_in)
%   % write_data(W_stable_in)
%   %
%   % Writes the stable manifold data to file.
% 
%   % Open file to write to
%   fileID = fopen('./data_files/2D_stable_manifold.txt', 'w');
%   % Write data
%   fprintf(fileID, '%.12E        %.12E        %.12E \n', [W_stable_in(:, 1), W_stable_in(:, 2), W_stable_in(:, 3)]);
% 
% end

function plot_stable_manifold(run_in, label_in)
  % plot_stable_manifold(run_in, label_in)
  %
  % Calculates, writes, and plots the 1-dimensional stable manifold of the
  % "central" saddle point, along with the periodic orbit.

  %-----------------------------------%
  %     Calculate Stable Manifold     %
  %-----------------------------------%
  % Base periodic orbit solution
  sol_PO = coll_read_solution('initial_PO', run_in, label_in);
  x_PO   = sol_PO.xbp;

  % Saddle-node point
  sol_EP = ep_read_solution('xpos', run_in, label_in);
  x_ep   = sol_EP.x;

  % Calculate stable manifold
  W_stable_p = calc_stable_manifold(run_in, label_in);

  % Write data to file
  % write_data(W_stable_p);

  %-------------------------------------------------------------------------%
  %%          Nonlinear Photonics 2024 Abstract - Periodic Orbit           %%
  %-------------------------------------------------------------------------%
  % Default line colours
  colours = colororder();

  % Setup figure
  fig = figure(1); clf;
  fig.Name = 'Initial Periodic Orbits (Phase)';
  fig.Units = 'inches'; fig.Position = [3, 3, 8, 7]; fig.PaperSize = [8, 7];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;
  ax.FontSize = 18;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plot original periodic orbit
  plot3(ax, x_PO(:, 1), x_PO(:, 2), x_PO(:, 3), Color=colours(3, :), ...
        LineWidth=4.0, DisplayName='$\Gamma$');

  % Plot equilibrium point
  plot3(ax, x_ep(1), x_ep(2), x_ep(3), Color='k', Marker='*', MarkerSize=20, ...
        MarkerFaceColor='k', MarkerEdgecolor='k', ...
        HandleVisibility='off');

  % Plot stable manifold
  plot3(ax, W_stable_p(:, 1), W_stable_p(:, 2), W_stable_p(:, 3), ...
        Color=colours(1, :), LineWidth=4.0, DisplayName='$W^{s}(p)$');

  % Hold axes
  hold(ax, 'off');

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 6.0];
  ax.YAxis.Limits = [0.0, 4.0];
  ax.ZAxis.Limits = [0.0, 21];

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickDirection = 'in';
  ax.XAxis.TickValues = 0.0 : 2.0 : 6.0;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 1.0 : 2.0 : 6.0;

  % Y-Axis
  ax.YAxis.TickDirection = 'in';
  ax.YAxis.TickValues = 0.0 : 1.0 : 4.0;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = 0.5 : 1.0 : 4.0;

  % Z-Axis
  ax.ZAxis.TickDirection = 'in';
  ax.ZAxis.TickValues = 0.0 : 5.0 : 20.0;
  ax.ZAxis.MinorTick = 'on';
  ax.ZAxis.MinorTickValues = 2.5 : 5.0 : 20.0;

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$G$';
  ax.YAxis.Label.String = '$Q$';
  ax.ZAxis.Label.String = '$I$';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % Grid lines
  ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

  % 3D plot view
  view(45, 10.0);

end
