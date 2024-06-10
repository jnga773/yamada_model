% Clear plots
close all;

% Clear workspace
clear;

% Add equation/functions to path
addpath('./functions/');
addpath('./plotting_scripts/');

%-------------------------------------------------------------------------%
%%             YAMADA MODEL (ode45 Phase-Space Trajectories)             %%
%-------------------------------------------------------------------------%
% Here we plot phase-space diagrams for trajectories of the Yamada model
% (using ode45) for some different chosen variables.
% The function containing the system of equations is written to
% './functions/yamada_ode.m'.

%--------------------%
%     PARAMETERS     %
%--------------------%
% Decay time of gain
gamma = 0.15;
% Pump current on the gain
A = 6.0;
% (Relative) absorption
B = 5.8;
a = 1.8;

% Time step
dt = 0.01;
% Max time to integrate over
t_max = 50.0;
% Set of times to integrate to steady state over
t = 0.0:dt:t_max;

% Initial conditions values
dx = 0.5;
x_max = 3.0;
x0_values = 0.0:dx:x_max;

% % Gain
% G0 = 0.0;
% % Absorption
% Q0 = 0.0;
% % Laser intensity
% I0 = 0.0;
% x0 = [G0; Q0; I0];

%-----------------------%
%     Solve (ode45)     %
%-----------------------%
% Empty arrays to add data
G_data = zeros(length(x0_values), length(x0_values), length(x0_values), length(t));
Q_data = zeros(length(x0_values), length(x0_values), length(x0_values), length(t));
I_data = zeros(length(x0_values), length(x0_values), length(x0_values), length(t));

% Redefine ode function with set parameters
% Temporary function for ode45
func_temp = @(t_in, x_in) ODE_yamada(t_in, x_in, [gamma; A; B; a]);

% Cycle through initial condition values
for j = 1:length(x0_values)
    for k = 1:length(x0_values)
        % Print progress to console
        fprintf('j = %d, k = %d, l = ... \n', j, k);

        for l = 1:length(x0_values)

            % Set initial conditions
            x0 = [x0_values(j); x0_values(k); x0_values(l)];

            % Solve using ode45
            [~, x_temp] = ode45(func_temp, t, x0);

            % Save data
            G_data(j, k, l, :) = x_temp(:, 1);
            Q_data(j, k, l, :) = x_temp(:, 2);
            I_data(j, k, l, :) = x_temp(:, 3);
        end
    end
end

%-------------------------------------------------------------------------%
%%                     PLOT PHASE SPACE TRAJECTORIES                     %%
%-------------------------------------------------------------------------%
% Run plotting script
plot_yamada_ode45_phase_space;
