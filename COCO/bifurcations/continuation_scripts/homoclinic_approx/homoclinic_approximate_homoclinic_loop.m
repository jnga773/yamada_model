%-------------------------------------------------------------------------%
%%              Reconstruct an Approximate Homoclinic Loop               %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_new = run_names.approx_homo.continue_homoclinics;
% Which run this continuation continues from
run_old = run_names.approx_homo.large_period_PO;

% Grab the label for the previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = max(label_old);

% Print to console
fprintf('~~~ Approximate Homoclinic: Third Run (ode_po2po) ~~~\n');
fprintf('Continue family of periodic orbits approximating homoclinics\n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s\n', label_old, run_old);

%--------------------------------------%
%     Initialise Problem Structure     %
%--------------------------------------%
% Initialise continuation problem
prob = coco_prob();

% Turn off bifurcation detections
prob = coco_set(prob, 'po', 'bifus', 'off');

% Set MXCL to false
prob = coco_set(prob, 'coll', 'MXCL', false);

% Continue a periodic orbit from a previous periodic orbit
prob = ode_po2po(prob, 'homo', run_old, label_old);

% Continue from equilibrium point
prob = ode_ep2ep(prob, 'x0', run_old, label_old);

% Glue parameters
prob = glue_parameters(prob);

% Assign 'gamma' to the set of active continuation parameters and 'po.period'
% to the set of inactive continuation parameters, thus ensuring that the
% latter is fixed during root finding.
prob = coco_xchg_pars(prob, 'gamma', 'homo.po.period');

% Set number of steps to confirm solution
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Set number of continuation steps
PtMX = 10000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Run continuation
coco(prob, run_new, [], 1, {'A', 'gamma', 'homo.po.period'}, p_range);

%-------------------------------------------------------------------------%
%%                               FUNCTIONS                               %%
%-------------------------------------------------------------------------%
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
