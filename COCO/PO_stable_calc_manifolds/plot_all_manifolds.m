%=========================================================================%
%                YAMADA MODEL (ATTRACTING PERIODIC ORBIT)                 %
%=========================================================================%
% In this script we plot the phase-space portrait of the attracting
% periodic orbit in the parameter region 7. Provided that the following
% scripts have been run,
%              - initial_periodic_orbit.m
%              - stable_manifold_q.m
%              - stable_manifold_PO.m

% Clear plots
close('all');

% Clear workspace
clear;
clc;

%-------------------------------------------------------------------------%
%%                               READ DATA                               %%  
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
% Load data from ./data_mat/initial_PO.mat
load('./data_mat/initial_PO.mat');

% Load stable manifold of q
load('./data_mat/stable_manifold_q.mat', 'Wsq');

% Load stable manifold of attracting periodic orbit
load('./data_mat/stable_manifold_PO.mat', 'WsPO_1', 'WsPO_2');

%-------------------%
%     Save Data     %
%-------------------%
data_out.xdim   = xdim;
data_out.pdim   = pdim;

data_out.p      = p;
data_out.pnames = pnames;

data_out.xbp_PO = xbp_PO;
data_out.tbp_PO = tbp_PO;
data_out.T_PO   = T;

data_out.x0     = x0;
data_out.xpos   = xpos;
data_out.xneg   = xneg;

data_out.Wsq    = Wsq;
data_out.WsPO_1 = WsPO_1;
data_out.WsPO_2 = WsPO_2;

% Save data to MATLAB .mat file
save('./data_mat/PO_and_manifolds.mat', '-struct', 'data_out');

%-------------------------------------------------------------------------%
%%                               PLOT DATA                               %%  
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

fig = figure(1); fig.Name = 'Initial Periodic Orbits'; clf;
fig.Units = 'inches'; fig.Position = [3, 3, 14, 8]; fig.PaperSize = [14, 8];

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

%------------------------------%
%     Plot: Periodic Orbit     %
%------------------------------%
% Hold axes
hold(ax, 'on');

% Plot stable periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
      LineStyle='-', Color=colours(3, :), ...
      DisplayName='$\Gamma$');

% Plot equilibrium points: x_{+}
plot3(ax, xpos(1), xpos(2), xpos(3), ...
      LineStyle='none', ...
      Marker='*', MarkerFaceColor='b', MarkerSize=10,  ...
      MarkerEdgeColor='b', DisplayName='$q$');

% Plot equilibrium points: x_{-}
plot3(ax, xneg(1), xneg(2), xneg(3), ...
      LineStyle='none', ...
      Marker='*', MarkerFaceColor='r', MarkerSize=10,  ...
      MarkerEdgeColor='r', DisplayName='$p$');

% Plot equilibrium points: x_{0}
plot3(ax, x0(1), x0(2), x0(3), ...
      LineStyle='none', ...
      Marker='o', MarkerFaceColor='r', MarkerSize=10, ...
      MarkerEdgeColor='r', DisplayName='$o$');

%------------------------------------%
%     Plot: Stable Manifold of q     %
%------------------------------------%
plot3(ax, Wsq(:, 1), Wsq(:, 2), Wsq(:, 3), Color=colours(1, :), ...
        DisplayName='$W^{s}(q)$');

%-------------------------------------%
%     Plot: Stable Manifold of PO     %
%-------------------------------------%
% Plot Surf plot (woah!)
  surf(ax, WsPO_1{1}, WsPO_1{2}, WsPO_1{3}, FaceColor=colours(1, :), FaceAlpha=0.25, ...
       MeshStyle='column', LineStyle='none', DisplayName='$W^{ss}(\Gamma)$');
  surf(ax, WsPO_2{1}, WsPO_2{2}, WsPO_2{3}, FaceColor=colours(1, :), FaceAlpha=0.25, ...
       MeshStyle='column', LineStyle='none', HandleVisibility='off');

%----------------%
%     Legend     %
%----------------%
% Legend
legend(ax, 'Interpreter', 'latex')

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-2.0, 8.0];
ax.YAxis.Limits = [-4.0, 8.0];
ax.ZAxis.Limits = [0.0, 35.0];

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$G(t)$';
ax.YAxis.Label.String = '$Q(t)$';
ax.ZAxis.Label.String = '$I(t)$';

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
