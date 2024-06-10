%-------------------------------------------------------------------------%
%%                          Setup Lin's Method                           %%
%-------------------------------------------------------------------------%
%-------------------------------------%
%     Parameters for Lin's Method     %
%-------------------------------------%
% Parameters from approximate homoclinic (run8)
% A     = 6.51653;
% gamma = 0.1;

% Parameters on left of MX point
% A     = 6.13692;
% gamma = 0.170562;

% Parameters on right of MX point
% A     = 6.68928;
% gamma = 0.0776709;

% A = 6.7622;
% gamma = 6.68078e-2;

% Parameters AT the MX point
% A     = 6.257098306041695;
% gamma = 0.141627582028018;

%-----------------------------------%
%     Start from run08 Solution     %
%-----------------------------------%
% Label from solution we start from
data_bcs.label_approx = 1;
% data_bcs.label_approx = 46;

% Read parameters
[data_bcs.p0, ~] = read_approximate_homoclinic_solution(data_bcs.label_approx);

%-----------------------------------------------%
%     Set Parameters and Equilibrium Points     %
%-----------------------------------------------%
% % Set new parameter vector
% data_bcs.p0 = [gamma; A; B; a];

% Calculate x0_pos and x0_neg
[data_bcs.x0_pos, data_bcs.x0_neg] = non_trivial_ss(data_bcs.p0);

% Calculate non-trivial steady states
[vu, vs1, vs2, eigvals] = unstable_stable_eigenvectors(data_bcs.x0_neg, data_bcs.p0);

% Save to structure
data_bcs.vu  = vu;
data_bcs.vs1 = vs1;
data_bcs.vs2 = vs2;
data_bcs.lu  = eigvals(1);
data_bcs.ls1 = eigvals(2);
data_bcs.ls2 = eigvals(3);

%----------------------------------%
%     Setup Lin's Method Stuff     %
%----------------------------------%
% Initial distances from the equilibria, along the tangent spaces of the
% unstable and stable manifolds, to the initial points along the corresponding
% trajectory segments.
eps1 = 0.009;
eps2 = 0.005;

% eps1 = 0.0009;
% eps2 = 0.0005;

% eps1 = -eps1;

% Angle for stable vector component
theta0 = 5.969;

% Lin epsilons vector
epsilon0 = [eps1; eps2; theta0];

%--------------------------------------------%
%     Boundary Conditions data Structure     %
%--------------------------------------------%
% Add parameter names to data_bcs
data_bcs.pnames = pnames;

% Normal vector to hyperplane \Sigma (just the y-axis at x=0.5)
data_bcs.normal = [0, 0, 1];
% Intersection point for hyperplane
data_bcs.pt0 = [1.25; 0.66; data_bcs.x0_pos(3)];

% Store the stable and unstable equilibria points into data_bcs
data_bcs.equilib_pt = data_bcs.x0_neg;

% Initial time
data_bcs.t0 = 0;

% Unstable Manifold: Initial point
data_bcs.x_init_u = data_bcs.equilib_pt' + eps1 * vu';
% Unstable Manifold: Final point
data_bcs.x_final_u = data_bcs.pt0;

% Stable Manifold: Initial point
data_bcs.x_init_s = data_bcs.equilib_pt' + eps2 * (cos(theta0) * vs1' + sin(theta0) * vs2');
% Stable Manifold: Final point
data_bcs.x_final_s = data_bcs.pt0;

% Plot "approximate" homoclinic with steady and unsteady eiegenvectors
% plot_homoclinic_manifolds(run7, p0_L, save_figure);
