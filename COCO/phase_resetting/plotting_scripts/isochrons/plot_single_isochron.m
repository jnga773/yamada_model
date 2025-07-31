function plot_single_isochron(run_in, labels_to_plot)
  % plot_single_isochron(run_in, label_plot)
  %
  % Plots a single isochron from the isochron_test.m run.
  arguments
    run_in;
    labels_to_plot char = '';
  end

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Bifurcation data matrix
  bd_read = coco_bd_read(run_in);

  % Read data
  iso1 = coco_bd_col(bd_read, 'iso1');
  iso2 = coco_bd_col(bd_read, 'iso2');
  iso3 = coco_bd_col(bd_read, 'iso3');

  % Get coordinates for 'FP' points
  if ~isempty(labels_to_plot)
    FP_idxs = coco_bd_idxs(bd_read, labels_to_plot);
  end

  % % Read unperturbed periodic orbit data
  % load('./data_mat/PO_and_manifolds.mat', 'xbp_PO', 'Ws_q', 'Ws_PO1', 'Ws_PO2');

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

  % Plot base orbit and manifolds
  plot_base_periodic_orbit(ax);

  hold(ax, 'on');

  % Plot single isochron
  plot3(ax, iso1, iso2, iso3, Color='k', LineStyle='-', ...
        DisplayName='Isochron');
  
  if ~isempty(labels_to_plot)
    plot3(ax, iso1(FP_idxs), iso2(FP_idxs), iso3(FP_idxs), ...
          MarkerFaceColor='r', LineStyle='none', Marker='o', ...
          DisplayName='FP');
  end

  hold(ax, 'off');

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
  % ax.XAxis.Limits = [-2.0, 6.0];
  % ax.YAxis.Limits = [-2.0, 6.0];
  % ax.ZAxis.Limits = [0.0, 50.0];

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