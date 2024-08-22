function plot_solutions_scan(run_in)
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Stable periodic orbit
  [sol_s, ~] = coll_read_solution('PO_stable', run_in, 1);
  xbp_PO_s = sol_s.xbp;

  % Unstable periodic orbit
  [sol_u, ~] = coll_read_solution('PO_unstable', run_in, 1);
  xbp_PO_u = sol_u.xbp;

  % Equilibrium points
  sol_0 = ep_read_solution('x0', run_in, 1);
  x_0   = sol_0.x;
  sol_pos = ep_read_solution('xpos', run_in, 1);
  x_pos = sol_pos.x;
  sol_neg = ep_read_solution('xneg', run_in, 1);
  x_neg = sol_neg.x;

  % Create empty data arrays
  XX1 = []; XX2 = [];
  YY1 = []; YY2 = [];
  ZZ1 = []; ZZ2 = [];

  % Cycle through stable manifold solutions
  bd = coco_bd_read(run_in);
  for label = coco_bd_labs(bd)
    % Grab solution
    [sol1, ~] = coll_read_solution('W1', run_in, label);
    [sol2, ~] = coll_read_solution('W2', run_in, label);
    
    % Append to data arrays
    XX1 = [XX1, sol1.xbp(:, 1)];
    YY1 = [YY1, sol1.xbp(:, 2)];
    ZZ1 = [ZZ1, sol1.xbp(:, 3)];

    XX2 = [XX2, sol2.xbp(:, 1)];
    YY2 = [YY2, sol2.xbp(:, 2)];
    ZZ2 = [ZZ2, sol2.xbp(:, 3)];
  end

  %-------------------%
  %     Plot Data     %
  %-------------------%
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

  % Plot stable periodic orbit
  plot3(ax, xbp_PO_s(:, 1), xbp_PO_s(:, 2), xbp_PO_s(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma_{s}$');

  % Plot unstable periodic orbit
  plot3(ax, xbp_PO_u(:, 1), xbp_PO_u(:, 2), xbp_PO_u(:, 3), ...
        LineStyle='-', Color=colours(4, :), ...
        DisplayName='$\Gamma_{u}$');

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

  % Plot Surf plot (woah!)
  surf(ax, XX1, YY1, ZZ1, FaceColor=colours(1, :), FaceAlpha=0.5, ...
       MeshStyle='column', LineStyle='-', EdgeColor=[0.6, 0.6, 0.6], ...
       LineWidth=0.5, DisplayName='$W^{s}(\Gamma_{\times})$');
  surf(ax, XX2, YY2, ZZ2, FaceColor=colours(1, :), FaceAlpha=0.5, ...
       MeshStyle='column', LineStyle='-', EdgeColor=[0.6, 0.6, 0.6], ...
       LineWidth=0.5, HandleVisibility='off');

  % Legend
  legend(ax, 'Interpreter', 'latex')

  hold(ax, 'off');

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 10.0];
  ax.YAxis.Limits = [0.0, 8.0];
  ax.ZAxis.Limits = [0.0, ceil(max(xbp_PO_s(:, 3)))];

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

end