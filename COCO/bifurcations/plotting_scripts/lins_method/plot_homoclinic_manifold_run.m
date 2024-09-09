function plot_homoclinic_manifold_run(run_in, label_in, data_in, fig_num_in)
  % PLOT_HOMOCLINIC_MANIFOLD_RUN: Plots the solution calculating in
  % run [run_in] for label [label_in].

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read solution of current run
  [sol1, ~] = coll_read_solution('unstable', run_in, label_in);
  [sol2, ~] = coll_read_solution('stable', run_in, label_in);

  % x-solution
  x_sol1 = sol1.xbp;
  x_sol2 = sol2.xbp;

  %-----------------------------------------------------------------------%
  %                    Plot: Homoclinic Manifold (3D)                     %
  %-----------------------------------------------------------------------%
  fig = figure(8); clf;
  fig.Name = 'Homoclinic Manifolds (3D)';
  fig.Units = 'inches';
  fig.Position = [2, 2, 12, 8];
  fig.PaperSize = [12, 8];

  % Setup axis
  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plot COCO solution
  plot3(ax, x_sol1(:, 1), x_sol1(:, 2), x_sol1(:, 3), LineStyle='-', ...
        Marker='.', MarkerSize=15, DisplayName='Unstable Manifold');
  plot3(ax, x_sol2(:, 1), x_sol2(:, 2), x_sol2(:, 3), LineStyle='-', ...
        Marker='.', MarkerSize=15, DisplayName='Stable Manifold');

  % Plot base solution
  plot_homoclinic_hyperplane_base(ax, run_in, label_in, data_in)

  % Legend
  legend(ax, 'Interpreter', 'latex')

  hold(ax, 'off');

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0, 7];
  ax.YAxis.Limits = [0, 7];
  ax.ZAxis.Limits = [0, 7];
  
  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$G(t)$';
  ax.YAxis.Label.String = '$Q(t)$';
  ax.ZAxis.Label.String = '$I(t)$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  title_str = sprintf('COCO Solution (run: $\\verb!%s!$, label: %d)', run_in, label_in);
  ax.Title.String = title_str;

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  % 3D view
  view(45, 15);

  %---------------------%
  %     Save Figure     %
  %---------------------%
  % if save_figure == true
  %   % Filename
  %   figname = sprintf('homoclinic_trajectory_%s', run_in(1:5));
  %   exportgraphics(fig, ['./images/', figname, '.pdf'], ContentType='vector');
  % end

end