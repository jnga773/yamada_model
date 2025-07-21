function data_out = calc_initial_solution_PR(run_in, label_in)
  % data_out = calc_initial_solution_PR(run_in, label_in)
  %
  % Reads data from previous run solution and calculates the 
  % initial conditions for the various different trajectory segments.
  %
  % Parameters
  % ----------
  % run_in : string
  %     The run identifier for the continuation problem.
  % label_in : integer
  %     The label identifier for the continuation problem.
  %
  % Returns
  % -------
  % data_out : struct
  %     Structure containing the initial conditions for the trajectory segments.
  %     Fields:
  %         - xdim : Original dimension of state space.
  %         - pdim : Original dimension of parameter space.
  %         - p0 : Initial parameter array.
  %         - pnames : Parameter names.
  %         - p_maps : Index mapping of each parameter.
  %         - t_seg1 : Initial time solutions for segment 1.
  %         - t_seg2 : Initial time solutions for segment 2.
  %         - t_seg3 : Initial time solutions for segment 3.
  %         - t_seg4 : Initial time solutions for segment 4.
  %         - x_seg1 : Initial state space solutions for segment 1.
  %         - x_seg2 : Initial state space solutions for segment 2.
  %         - x_seg3 : Initial state space solutions for segment 3.
  %         - x_seg4 : Initial state space solutions for segment 4.
  %
  % See Also
  % --------
  % coll_read_solution

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    label_in double
  end

  %-----------------------------------------------------------------------%
  %                            Read Data                                  %
  %-----------------------------------------------------------------------%
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read COCO solution
  [sol, data] = coll_read_solution('initial_PO', run_in, label_in);

  % Original dimension of state space
  xdim = data.xdim;
  % Original dimension of parameter space
  pdim = data.pdim;

  % State space solution
  xbp_read = sol.xbp;
  
  % Periodic orbit solution
  gamma_read = xbp_read(:, 1:xdim);

  % Initial zero-phase point of the periodic orbit
  gamma_0 = gamma_read(1, :)';

  % Time data
  tbp_read = sol.tbp;

  % Read parameters
  p_read      = sol.p;
  pnames_read = data.pnames;

  % Calculate stationary point
  [xpos, ~] = non_trivial_ss(p_read(1:pdim));

  %--------------------%
  %     Parameters     %
  %--------------------%
  % System parameters
  p_system      = p_read(1:pdim);
  pnames_system = {};
  for i = 1 : pdim
    pnames_system{i} = pnames_read{i};
  end

  %----------------------------%
  %     Initial Parameters     %
  %----------------------------%
  % Phase along the periodic orbit
  theta_gamma = 1.0;

  %---------------------------%
  %     Parameter Indices     %
  %---------------------------%
  % This p_maps data structure will be used in each boundary condition
  % function to ensure the correct parameters are being used.

  % Save the index mapping of each parameter
  p_maps.theta_gamma = pdim + 1;

  %------------------------%
  %     Set Parameters     %
  %------------------------%
  % Initial parameter array
  p0_out = zeros(pdim+length(p_maps), 1);
  % Put parameters in order
  p0_out(1:pdim)             = p_system;
  p0_out(p_maps.theta_gamma) = theta_gamma;

  %-------------------------%
  %     Parameter Names     %
  %-------------------------%
  % Parameter names
  pnames_PR                     = {pnames_system{1:pdim}};
  pnames_PR{p_maps.theta_gamma} = 'theta_gamma';

  %-----------------------------------------------%
  %     Segment Initial Conditions: Easy Mode     %
  %-----------------------------------------------%
  % Segment 1
  t_seg1 = tbp_read;
  x_seg1 = gamma_read;
  
  % Segment 2
  t_seg2 = [0.0; max(tbp_read)];
  x_seg2 = [gamma_0'; gamma_0'];

  %----------------%
  %     Output     %
  %----------------%
  % Original vector field dimensions
  data_out.xdim       = xdim;
  data_out.pdim       = pdim;

  % Parameters
  data_out.p0         = p0_out;
  data_out.pnames     = pnames_PR;
  data_out.p_maps     = p_maps;

  % Initial time solutions for each segment
  data_out.t_seg1     = t_seg1;
  data_out.t_seg2     = t_seg2;

  % Initial state space solutions for each segment
  data_out.x_seg1     = x_seg1;
  data_out.x_seg2     = x_seg2;

  % Equilibrium point
  data_out.xpos       = xpos;

end