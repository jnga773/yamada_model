function plot_orbit_and_Wq_solution(run_in, label_in)
  % plot_orbit_and_Wq_solution(run_in, label_in)
  %
  % Plots the trajectory segment of the stable manifold of q

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read phase resetting solutions
  [sol1, ~] = coll_read_solution('seg1', run_in, label_in);
  [sol2, ~] = coll_read_solution('seg2', run_in, label_in);

  % Get state space solutions
  xbp1 = sol1.xbp;
  xbp2 = sol2.xbp;

  % Equilibrium point
  sol_pos = ep_read_solution('xpos', run_in, label_in);
  xpos = sol_pos.x;

  % Read manifold trajectory segment solution
  [solW, ~] = coll_read_solution('Wsq', run_in, label_in);
  % Get state space solutions
  xbpW = solW.xbp;

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

  % Plot phase resetting segment 1
  plot3(ax, xbp1(:, 1), xbp1(:, 2), xbp1(:, 3), ...
        LineStyle='-', Color=colours(2, :), ...
        DisplayName='seg1');

  % Plot phase resetting segment 2
  plot3(ax, xbp2(:, 1), xbp2(:, 2), xbp2(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='seg2');

  % Plot equilibrium points: x_{+}
  plot3(ax, xpos(1), xpos(2), xpos(3), ...
        LineStyle='none', ...
        Marker='o', MarkerFaceColor='r', MarkerSize=10,  ...
        MarkerEdgeColor='r', DisplayName='$q$');

  % Plot manifold trajectory segment
  plot3(ax, xbpW(:, 1), xbpW(:, 2), xbpW(:, 3), Color=colours(1, :), ...
        DisplayName='$W^{s}(q)$');

  legend(ax, Interpreter='latex');

  hold(ax, 'off');

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 5.0];
  ax.YAxis.Limits = [0.0, 4.0];
  ax.ZAxis.Limits = [0.0, 16];

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