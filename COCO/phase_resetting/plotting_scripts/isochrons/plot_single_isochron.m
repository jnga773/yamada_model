function plot_single_isochron(run_in)
  % plot_single_isochron(run_in, save_figure)
  %
  % Plots a single isochron from the isochron_test.m run.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Bifurcation data matrix
  bd_read = coco_bd_read(run_in);

  % Read data
  iso1 = coco_bd_col(bd_read, 'iso1');
  iso2 = coco_bd_col(bd_read, 'iso2');
  iso3 = coco_bd_col(bd_read, 'iso3');

  % Read unperturbed periodic orbit data
  load('./data_mat/initial_PO.mat');

  %-------------------------------------------------------------------------%
  %%                         Plot: Single Isochron                         %%
  %-------------------------------------------------------------------------%
  % matplotlib colour order
  colours = colororder();

  fig = figure(1);
  fig.Name = 'Single (Test) Isochron';
  fig.Units = 'inches';
  fig.Position = [3, 3, 8, 8]; fig.PaperSize = [8, 8];

  tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
  ax = nexttile;

  %--------------%
  %     Plot     %
  %--------------%
  hold(ax, 'on');

  % Plot single isochron
  plot3(ax, iso1, iso2, iso3, Color='k', LineStyle='-', ...
        DisplayName='Isochron');

  % Plot unperturbed periodic orbit
  plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
        LineStyle='-', Color=colours(3, :), ...
        DisplayName='$\Gamma$');

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
  % % X-Axis
  % ax.XAxis.TickValues = 0.0 : 1.0 : 5.0;
  % ax.XAxis.MinorTick = 'on';
  % ax.XAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

  % % Y-Axis
  % ax.YAxis.TickValues = 0.0 : 1.0 : 5.0;
  % ax.YAxis.MinorTick = 'on';
  % ax.YAxis.MinorTickValues = 0.5 : 1.0 : 5.0;

  % % Z-Axis
  % ax.ZAxis.TickDirection = 'in';
  % % ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
  % % ax.ZAxis.MinorTick = 'on';
  % % ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  % ax.XAxis.Limits = [0.0, 5.0];
  % ax.YAxis.Limits = [0.0, 4.0];
  % ax.ZAxis.Limits = [0.0, ceil(max(xbp_PO(:, 3)))];

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

end