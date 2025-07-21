function plot_initial_PO(run_in, label_in)
  % plot_initial_PO(run_in, label_in)
  %
  % Plots the initial periodic orbit from the 'coll' toolbox run.
  % along with the three stationary points ('q', 'p', and 'o'),
  % and the one-dimensional stable manifold of point 'q'.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read periodic orbit solution
  [sol_PO, data_PO] = coll_read_solution('initial_PO', run_in, label_in);
  % Periodic orbit
  xbp_PO = sol_PO.xbp;
  % tbp_PO = sol_PO.tbp;
  % T_PO   = sol_PO.T;

  % Read equilibrium points
  sol_0   = ep_read_solution('x0', run_in, label_in);
  sol_neg = ep_read_solution('xneg', run_in, label_in);
  sol_pos = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium points
  x0   = sol_0.x;
  xpos = sol_pos.x;
  xneg = sol_neg.x;

  % Calculate stable manifold of x_{+} (q)
  Wq_s = calc_stable_manifold(run_in, label_in);

  %--------------------------------------%
  %     Plot Initial Periodic Orbits     %
  %--------------------------------------%
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

  % Plot initial periodic orbit
  plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
        LineStyle='-', LineWidth=2.0, Color=colours(3, :), ...
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
  
  % Plot stable manifold of q / x_{+}
  plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
        Color=colours(1, :), LineWidth=2.0, ...
        DisplayName='$W^{s}(p)$');

  % Legend
  legend(ax, Interpreter='latex');

  % Turn of axis hold
  hold(ax, 'off');

  %--------------------%
  %     Axis Ticks     %
  %--------------------%
  % X-Axis
  ax.XAxis.TickDirection = 'in';
  ax.XAxis.TickValues = 0.0 : 2.0 : 10.0;
  ax.XAxis.MinorTick = 'on';
  ax.XAxis.MinorTickValues = 1.0 : 2.0 : 10.0;
  
  % Y-Axis
  ax.YAxis.TickDirection = 'in';
  ax.YAxis.TickValues = 0.0 : 2.0 : 8.0;
  ax.YAxis.MinorTick = 'on';
  ax.YAxis.MinorTickValues = 1.0 : 2.0 : 8.0;
  
  % Z-Axis
  ax.ZAxis.TickDirection = 'in';
  % ax.ZAxis.TickValues = 0.0 : 2.0 : 18.0;
  % ax.ZAxis.MinorTick = 'on';
  % ax.ZAxis.MinorTickValues = 1.0 : 2.0 : 18.0;

  %---------------------%
  %     Axis Limits     %
  %---------------------%
  ax.XAxis.Limits = [0.0, 10.0];
  ax.YAxis.Limits = [0.0, 8.0];
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

end

function x_out = calc_stable_manifold(run_in, label_in)
      % x_out = calc_stable_manifold(run_in, label_in)
      %
      % Calculate the stable manifold of the equilibrium point in the middle of
      % the periodic orbit.
    
      %-------------------%
      %     Arguments     %
      %-------------------%
      arguments
        run_in char
        label_in double
      end
    
      %-------------------%
      %     Read Data     %
      %-------------------%
      % Read EP solution
      sol_EP   = ep_read_solution('xpos', run_in, label_in);
    
      % Equilibrium point
      x_pos = sol_EP.x;
      % Parameters
      p     = sol_EP.p;
      
      %------------------------------%
      %     Calculate EigenStuff     %
      %------------------------------%
      % Jacobian
      J_stable = yamada_DFDX(x_pos, p);
    
      % Calculate eigenvalues and eigenvectors
      [eigvec, eigval] = eig(J_stable);
    
      % Indices for stable eigenvectors (eigval < 0)
      % stable_index = find(diag(eigval) < 0);
      stable_index = 3;
    
      % Stable eigenvector
      vec_s = eigvec(:, stable_index);
    
      %------------------------------------------------%
      %     Calculate Things: Positive y direction     %
      %------------------------------------------------%
      % Small distance
      eps1 = -0.01;
      % Time span
      % t_span1 = [0.0, -17.0];
      t_span1 = -16.5 : 0.01 : 0.0;
      t_span1 = flip(t_span1);
    
      % Initial vector
      x_init1 = x_pos + (eps1 * vec_s);
    
      % Integrate using ode45
      [~, W1] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span1, x_init1);
    
      %------------------------------------------------%
      %     Calculate Things: Positive y direction     %
      %------------------------------------------------%
      % Small distance
      eps2 = 0.01;
      % Time span
      % t_span2 = [0.0, -25.0];
      t_span2 = -25.0 : 0.01 : 0.0;
      t_span2 = flip(t_span2);
    
      % Initial vector
      x_init2 = x_pos + (eps2 * vec_s);
    
      % Integrate using ode45
      [~, W2] = ode45(@(t_in, x0_in) yamada(x0_in, p), t_span2, x_init2);
    
      %----------------%
      %     Output     %
      %----------------%
      x_out = [flip(W2); W1];
    
    end