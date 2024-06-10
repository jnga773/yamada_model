function plot_homoclinic_manifolds(p0_in, save_figure)
  % PLOT_HOMOCLINIC_MANIFOLDS: Plots the solution calculating in run [run_in] for label
  % [label_in].

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Calculate non-trivial steady states and stable and unstable eigenvectors
  [x0_neg, vu, vs1, vs2] = unstable_stable_eigevectors(p0_in);

  % Scale for vectors
  scale = 0.5;

  % Shift vectors by the stationary point
  vu = [x0_neg, x0_neg + scale * vu];
  vs1 = [x0_neg, x0_neg - scale * vs1];
  vs2 = [x0_neg, x0_neg + scale * vs2];

  %-----------------------------------------------------------------------%
  %                    Plot: Homoclinic Manifold (3D)                     %
  %-----------------------------------------------------------------------%
  fig = figure(13); clf;
  fig.Name = 'Homoclinic Manifolds (3D)';
  fig.Units = 'inches';
  fig.Position = [0, 0, 12, 8];
  fig.PaperSize = [12, 8];

  % Setup axis
  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  % Hold axes
  hold(ax, 'on');

  % Plot base solution
  plot_homoclinic_hyperplane_base(ax, p0_in)

  % Legend
  legend(ax, 'Interpreter', 'latex')

  hold(ax, 'off');

  %---------------------%
  %     Axis Labels     %
  %---------------------%
  ax.XAxis.Label.String = '$G(t)$';
  ax.YAxis.Label.String = '$Q(t)$';
  ax.ZAxis.Label.String = '$I(t)$';

  %--------------------%
  %     Axis Title     %
  %--------------------%
  ax.Title.String = 'Homoclinic Manifolds';

  %----------------------%
  %     Figure Stuff     %
  %----------------------%
  box(ax, 'on');
  grid(ax, 'on');

  %---------------------%
  %     Save Figure     %
  %---------------------%
  if save_figure == true
    exportgraphics(fig, './images/homoclinic_manifolds.pdf', ContentType='vector');
  end

end