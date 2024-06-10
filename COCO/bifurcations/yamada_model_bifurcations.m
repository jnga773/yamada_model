% Clear plots
close('all');

% Clear workspace
clear;
clc;

% Add equation/functions to path
addpath('./functions/');
addpath('./continuation_scripts/');
addpath('./plotting_scripts/');

% Save figures switch
% save_figure = true;
save_figure = false;

%-------------------------------------------------------------------------%
%%               YAMADA MODEL (COCO Bifurcation Diagrams)                %%
%-------------------------------------------------------------------------%
% We compute a family of bifurcations from a driven-pulsed laser, described
% by the Yamada model of coupled differential equations:
%     G' = \gamma (A - G - G I) ,
%     Q' = \gamma (B - Q - a Q I) ,
%     I' = (G - Q - 1) I ,
% where G is the gain, Q is the absorption, and I is the intensity of the
% laser. The system is dependent on four parameters: the pump current on
% the gain, A (or A); the relative absoprtion, B and a; and the decay
% time of the gain, \gamma.

% Here we plot bifurcation diagrams of the Yamada model using COCO. The
% functions used are written t
%     - './functions/yamada.m' ,
%     - './functions/yamada_DFDX.m' , and
%     - './functions/yamada_DFDP.m' ,
% or, in the COCO Symbolic form:
%     - './functions/yamada_symbolic.m'.

%--------------------%
%     Parameters     %
%--------------------%
% Because we will only be looking at the (A, \gamma) plane, we will be
% setting values for a and B.
B = 5.8;
a = 1.8;

% Set some initial values for \gamma and A
gamma = 0.10;
A = 6.6;

%--------------------%
%     COCO Setup     %
%--------------------%
% Parameter names
pnames = {'gamma', 'A', 'B', 'a'};

% Initial parameter values
p0 = [gamma; A; B; a];

% Initial state values
x0 = [A; B; 0];

% Parameter ranges
gamma_range = [0.0, 0.25];
A_range = [5.0, 11.0];
p_range = {A_range, gamma_range};

% State dimensions
pdim = length(p0);
xdim = length(x0);

%------------------------------%
%     Analytic Expressions     %
%------------------------------%
% % Print some analytic state things
% fprintf('~~~ Some Analytical Results ~~~\n');
% 
% % Line of transcritical bifurcations
% A_T = B + 1;
% fprintf('Line of transcritical bifurcations at A_T = %.3f \n', A_T);
% 
% % Line of Saddle-Node bifurcations
% A_S = (1/a) * (-1 + a + B + 2 * sqrt(B * (a - 1)));
% fprintf('Line of saddle-node bifurcations at A_S = %.3f \n', A_S);
% 
% % Bogdanov-Takens points
% A_BT = A_S;
% gamma_BT = -(1 + (B * (a - 1)) - 2 * sqrt(B * (a - 1))) / (sqrt(B * (a - 1)) * (1 - a - sqrt(B * (a - 1))));
% fprintf('Bogdanov-Takens points at (A_BT, g_BT) = (%.3f, %.3f) \n', A_BT, gamma_BT);
% 
% fprintf('\n');

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Hard-coded function list
% funcs = {@yamada, @yamada_DFDX, @yamada_DFDP};
funcs = yamada_symbolic();

%-------------------------------------------------------------------------%
%%                         Initial Continuation                          %%
%-------------------------------------------------------------------------%
% We set up the continuation problem by following the equilibrium point x0.
% We then follow the branching point, from which many of the other 
% bifurcations come about.

%-------------------------------------%
%%     Compute Equilibrium Point     %%
%-------------------------------------%
% Find initial equilibrium points along line \gamma = 0.08
run_names.initial_run = 'run01_initial_run';

% Run continuation
% initial_continuation;

%------------------------------------%
%%     Continue Branching Point     %%
%------------------------------------%
% Continue from the branching point found in run1
run_names.branching_point = 'run02_branching_point_continuation';

% Run continuation
% equilibrium_from_branching_point;

%------------------------------------%
%%     Hopf Bifurcation Line: H     %%
%------------------------------------%
% Continue the Hopf bifurcation point
run_names.hopf_bifurcations = 'run03_hopf_bifurcation_line_H';

% Run continuation
% hopf_bifurcation_line;

%---------------------------------%
%%     Saddle-Node Line: A_S     %%
%---------------------------------%
% Continue the saddle-node point from run2 with a two parameter
% continuation to find the saddle-node line at A_S.
run_names.saddle_nodes = 'run04_saddle_node_line_AS';

% Run continuation
% saddle_node_line;

%-----------------------------------%
%%     Transcritical Line: A_T     %%
%-----------------------------------%
% Continue the equilibrium point from run2 at the A_T line.
run_names.transcritical = 'run05_transcritical_points_line_AT';

% Run continuation
% transcritical_line;

%-------------------------------------------------------------------------%
%%                Compute Homoclinic Orbits (Approximate)                %%
%-------------------------------------------------------------------------%
% We compute a family of homoclinic orbits emanating from a Hopf
% bifurcation. We approximate the homoclinic orbit as a periodic orbit with
% very large period.

% Add approximate homoclinic scripts to path
addpath('./continuation_scripts/homoclinic_approx/');
addpath('./plotting_scripts/homoclinic_approx/');

%---------------------------------------------%
%%     Compute Family of Periodic Orbits     %%
%---------------------------------------------%
% Continue from Hopf bifurcation to periodic orbit
run_names.approx_homo.PO_from_hopf = 'run06_periodic_orbit_continuation_from_hopf_bifurcation';

% Run continuation
% homoclinic_hopf_to_periodic_orbit;

%----------------------------------%
%%     Large T Periodic Orbit     %%
%----------------------------------%
% Extend the periodic of periodic orbit from run06 and reconstruct the
% periodic orbit.
run_names.approx_homo.large_period_PO = 'run07_reconstruct_large_period_periodic_orbit';

% Run continuation
% homoclinic_reconstruct_high_period_periodic_orbit;

%------------------------------------------%
%%     Follow Family of "Homoclinics"     %%
%------------------------------------------%
% Follow family of "homoclinics" in A and \gamma.
run_names.approx_homo.continue_homoclinics = 'run08_approximate_homoclinic';

% Run continuation
% homoclinic_approximate_homoclinic_loop;

%-------------------------------------------------------------------------%
%%                      Compute Double Limit Cycle                       %%
%-------------------------------------------------------------------------%
% Compute and follow the double limit cycles (saddle-nodes of periodic
% orbits). 

% Add double limit cycle scripts to path
addpath('./continuation_scripts/double_limit_cycle/');
addpath('./plotting_scripts/double_limit_cycle/');

%--------------------------------------------%
%%     Calculate Initial Periodic Orbit     %%
%--------------------------------------------%
% Find initial periodic orbit
run_names.limit_cycle.initial_PO = 'run09_initial_periodic_orbit';

% Run continuation
% double_limit_cycle_initial_periodic_orbit;

%--------------------------------------%
%%     Follow Double Limit Cycles     %%
%--------------------------------------%
% Continue from periodic orbit to find folding point
run_names.limit_cycle.follow_limit_cycle = 'run10_double_limit_cycle';

% Run continuation
% double_limit_cycle_continue_from_saddle;

%-------------------------------------------------------------------------%
%%                Compute Homoclinic Orbits (Lin's Method)               %%
%-------------------------------------------------------------------------%
% We use Lin's method to calculate the homoclinic orbits. This is a better
% approximation than the large period periodic orbits, as above.

% Add Lin's method scripts to path
% Add script folders to path
addpath('./continuation_scripts/lins_method/');
addpath('./continuation_scripts/lins_method/boundary_conditions/');
addpath('./plotting_scripts/lins_method/');

%-------------------------------------%
%%     Compute Unstable Manifold     %%
%-------------------------------------%
% Grow unstable manifold
run_names.lins_method.unstable_manifold = 'run11_unstable_manifold';

% Run continuation
lins_method_unstable_manifold;

%-----------------------------------%
%%     Compute Stable Manifold     %%
%-----------------------------------%
% Grow stable manifold
run_names.lins_method.stable_manifold = 'run12_stable_manifold';

% Run continuation
lins_method_stable_manifold;

%-----------------------------%
%%     Close the Lin Gap     %%
%-----------------------------%
% Close Lin gap
run_names.lins_method.close_lingap = 'run13_close_lingap';

% Run continuation
lins_method_close_lingap;

%-------------------------------%
%%     Close Distance eps1     %%
%-------------------------------%
% Close distance eps1
run_names.lins_method.close_eps1 = 'run14_close_eps1';

% Run continuation
lins_method_close_eps1;

%-------------------------------%
%%     Close Distance eps2     %%
%-------------------------------%
% Close distance eps2
run_names.lins_method.close_eps2 = 'run15_close_eps2';

% Run continuation
lins_method_close_eps2;

%-----------------------------------------------%
%%     Compute Family of Homoclinic Orbits     %%
%-----------------------------------------------%
% Continue along family of homoclinics
run_names.lins_method.continue_homoclinics = 'run16_continue_homoclinics';

% Run continuation
lins_method_continue_homoclinics;
% % lins_method_single_segment;

%-------------------------------------------------------------------------%
%%                                PLOTS                                  %%
%-------------------------------------------------------------------------%
% Write data to txt files
write_bifurcation_data(run_names);

% Run the plotting file
plot_bifurcation_diagram(p0, run_names, save_figure);
plot_bifurcation_diagram_zoomed(p0, run_names, save_figure); 

% % Plot homoclinic method comparison diagram
% compare_homoclinic_bifurcations(run_names, save_figure);
% test_compare_homoclinics(run_names, true);
