%-------------------------------------------------------------------------%
%%         Locate a High-Period Periodic Orbit from Previous Run         %%
%-------------------------------------------------------------------------%
% We approximate a homoclinic orbit by finding a periodic orbit that has a
% large period. The closer you are to an equilibrium point, the longer you
% spend then, and hence you have a super duper large period.

% We obtain an initial solution guess by extending the duration spent near
% an equilibrium point for a periodic orbit found in the previous run.
% The call to coco_xchg_pars constrains the interval length, while
% releasing the second problem parameter.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.approx_homo.large_period_PO;
% Which run this continuation continues from
run_old = run_names.approx_homo.PO_from_hopf;

% Read solution of previous run with largest period.
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = max(label_old);

% Print to console
fprintf('~~~ Approximate Homoclinic: Second Run (ode_isol2po) ~~~\n');
fprintf('Find reconstructed high-period periodic orbit approximating a homoclinic connection\n');
fprintf('Run name: %s\n', run_new);
fprintf('Continuing from point %d in run: %s\n', label_old, run_old);

%-------------------------------------%
%     Read Data from Previous Run     %
%-------------------------------------%
% Find minimum point corresponding to equilibrium point, and insert
% large time segment.
po_data = insert_large_time_segment(run_old, label_old);

%------------------------------------------%
%     Continuation from Periodic Orbit     %
%------------------------------------------%
% Continuation with a zero-dimensional atlas algorithm can now be
% used to locate a periodic orbit with period equal to its initial
% value, i.e., to [scale] times the largest period found in the previous
% run.

% Initialize continuation problem structure with the same number of
% intervals as in previous run.
prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', po_data.NTST);
prob = coco_set(prob, 'po', 'bifus', 'off');

% The value of 10 for 'NAdapt' implied that the trajectory discretisation
% is changed adaptively ten times before the solution is accepted.
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Continue periodic orbit from initial solution
prob = ode_isol2po(prob, 'homo', funcs{:}, ...
                   po_data.t_sol, po_data.x_sol, pnames, po_data.p0);

% Add instance of equilibrium point to find and follow the actual 
% equilibrium point
prob = ode_isol2ep(prob, 'x0', funcs{:}, ...
                   po_data.x0, po_data.p0);

% Glue parameters
prob = glue_parameters(prob);

% Assign 'gamma' to the set of active continuation parameters and 'po.period'
% to the set of inactive continuation parameters, thus ensuring that the
% latter is fixed during root finding.
prob = coco_xchg_pars(prob, 'gamma', 'homo.po.period');

% bd7 = coco(prob, run_new, [], 0, {'A', 'po.orb.coll.err_TF', 'po.period'});
coco(prob, run_new, [], 0, {'A', 'homo.po.orb.coll.err_TF', 'homo.po.period'});

%-------------------------------------------------------------------------%
%%                               Test Plot                               %%
%-------------------------------------------------------------------------%
%---------------------------%
%     Quick Test Things     %
%---------------------------%
% Quick test plot
% plot_2D_increasing_period;
% plot_increasing_period_3D;
% plot_high_period_time;

% Print some things
% fprintf('Large Period = %f \n', T);
% fprintf('G(end) = %f \n', sol.xbp(end, 1));
% fprintf('Q(end) = %f \n', sol.xbp(end, 2));
% fprintf('I(end) = %f \n', sol.xbp(end, 3));

% Quick test plot
plot_extended_time_profile(run_new, 2, save_figure)
plot_new_solution(run_old, run_new, save_figure);

%-------------------------------------------------------------------------%
%%                               FUNCTIONS                               %%
%-------------------------------------------------------------------------%
function po_data_out = insert_large_time_segment(run_in, label_in)
  % po_data_out = insert_large_time_segment(run_in, label_in)
  %
  % Reads the periodic orbit solution from solution [label_old] of
  % [run_old], and finds the segment of the state-space solution
  % closest to the equilibrium point.
  %
  % With this point found, we insert a large time segment to
  % "trick" the periodic orbit into having a larger period.
  %
  % Input
  % ----------
  % run_in: string
  %     The string identifier for the previous COCO run that we will
  %     read information from.
  % label_in: int
  %     The solution label we will read the data from
  %
  % Output
  % ----------
  % po_data_out : data structure
  %     Contains the state space solution, temporal solution and
  %     parameters.

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read solution with maximum period
  [sol, data] = coll_read_solution('po.orb', run_in, label_in);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Evaluate vector field at basepoints
  f = yamada(sol.xbp', repmat(sol.p, [1, size(sol.xbp, 1)]));

  % Extract the discretisation points corresponding to the minimum value of
  % the norm of the vector field along the longest-period periodic orbit.
  % Find basepoint closest to equilibrium
  f_norm = sqrt(sum(f .* f, 1)); f_norm = f_norm';
  [~, idx] = min(f_norm);

  % Print
  % fprintf('\n');
  % fprintf('idx = %d \n', idx);

  % Then insert a time segment that is a large multiple of the orbit
  % period immediately following the discretisation point.
  scale = 1000;
  T = sol.T;

  % fprintf('Maximum Period from run ''%s'', T = %f \n', run_new, T);
  % fprintf('Scaled period is T'' = %d x %f = %f \n', scale, T, scale * T);

  % Crank up period by factor scale
  t_sol = [sol.tbp(1:idx,1);
           T * (scale - 1) + sol.tbp(idx+1:end,1)];

  % Approximate equilibrium point
  x0 = sol.xbp(idx, :);

  %----------------%
  %     Output     %
  %----------------%
  % Temporal solution
  po_data_out.t_sol = t_sol;
  % State-space solution
  po_data_out.x_sol = sol.xbp;
  % Parameters
  po_data_out.p0    = sol.p;
  % Equilibrium point
  po_data_out.x0    = x0;
  % NTST setting from previous run
  po_data_out.NTST  = data.coll.NTST;

end

function prob_out = glue_parameters(prob_in)
  % prob_out = glue_parameter(prob_in)
  %
  % Glue the parameters of the EP segments and PO segment together 
  % (as they're all the same anyway)

  %---------------%
  %     Input     %
  %---------------%
  % Input continuation problem structure
  prob = prob_in;

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read index data periodic orbit segment
  [data, uidx] = coco_get_func_data(prob, 'homo.po.orb.coll', 'data', 'uidx');

  % Read index data equilibrium points
  [data1, uidx1] = coco_get_func_data(prob, 'x0.ep', 'data', 'uidx');

  % Index mapping
  maps = data.coll_seg.maps;
  maps1 = data1.ep_eqn;

  %-------------------------%
  %     Glue Parameters     %
  %-------------------------%
  prob = coco_add_glue(prob, 'shared_parameters', ...
                      uidx(maps.p_idx), uidx1(maps1.p_idx));

  %----------------%
  %     Output     %
  %----------------%
  prob_out = prob;

end
