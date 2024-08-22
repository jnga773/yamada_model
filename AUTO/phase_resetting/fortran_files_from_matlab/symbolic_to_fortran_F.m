% State space dimension
xdim = 3;

%-------------------------------------------------------------------------%
%%                      Phase Resetting Components                       %%
%-------------------------------------------------------------------------%
%-------------------%
%     Segment 1     %
%-------------------%
% State vectors
x_seg1 = sym('x1_vec', [xdim, 1]);
% Adjoint vetors
w_seg1 = sym('w1_vec', [xdim, 1]);

%-------------------%
%     Segment 2     %
%-------------------%
% State vectors
x_seg2 = sym('x2_vec', [xdim, 1]);
% Adjoint vetors
w_seg2 = sym('w2_vec', [xdim, 1]);

%-------------------%
%     Segment 3     %
%-------------------%
% State vectors
x_seg3 = sym('x3_vec', [xdim, 1]);

%-------------------%
%     Segment 4     %
%-------------------%
% State vectors
x_seg4 = sym('x4_vec', [xdim, 1]);

%--------------------%
%     Parameters     %
%--------------------%
% Yamada model parameters
syms gam A_pump B a
p_sys = [gam; A_pump; B; a];

% Phase resetting parameters
syms T k
syms theta_old theta_new
syms mu_s eta
syms A_perturb theta_perturb phi_perturb
p_PR = [T; k; theta_old; theta_new;
        mu_s; eta;
        A_perturb; theta_perturb; phi_perturb];

%-------------------------------------------------------------------------%
%%                         Vector Field Encoding                         %%
%-------------------------------------------------------------------------%
%-------------------%
%     Segment 1     %
%-------------------%
% Vector field
vec_field = yamada(x_seg1, p_sys);

% Vector equation
vec_eqn = T * theta_new * vec_field;

% Calculate tranpose of Jacobian at point xvec
J_T = transpose(jacobian(vec_field, x_seg1));

% Adjoint equation
adj_eqn = -T * theta_new * J_T * w_seg1;

% Total equation
F_seg1 = [vec_eqn; adj_eqn];

%-------------------%
%     Segment 2     %
%-------------------%
% Vector field
vec_field = yamada(x_seg2, p_sys);

% Vector equation
vec_eqn = T * (1 - theta_new) * vec_field;

% Calculate tranpose of Jacobian at point xvec
J_T = transpose(jacobian(vec_field, x_seg2));

% Adjoint equation
adj_eqn = -T * (1 - theta_new) * J_T * w_seg2;

% Total equation
F_seg2 = [vec_eqn; adj_eqn];

%-------------------%
%     Segment 3     %
%-------------------%
% Vector field
vec_field = yamada(x_seg3, p_sys);

% Vector equation
vec_eqn = T * (1 - theta_old) * vec_field;

% Total equation
F_seg3 = vec_eqn;

%-------------------%
%     Segment 4     %
%-------------------%
% Vector field
vec_field = yamada(x_seg4, p_sys);

% Vector equation
vec_eqn = T * k * vec_field;

% Total equation
F_seg4 = vec_eqn;

%-----------------------%
%     Write to File     %
%-----------------------%
% Full vector field
F_eqn = [F_seg1; F_seg2; F_seg3; F_seg4];

% Write to file
fortran(F_eqn, File='./fortran_generated/Segs_F.f90');

%-------------------------------------------------------------------------%
%%                               JACOBIANS                               %%
%-------------------------------------------------------------------------%
% Full state vector
uvec = [x_seg1; w_seg1; x_seg2; w_seg2; x_seg3; x_seg4];

% Full parameter vector
pvec = [p_sys; p_PR];

%------------------------------%
%     State-Space Jacobian     %
%------------------------------%
DFDX = jacobian(F_eqn, uvec);

% Write to file
fortran(DFDX, File='./fortran_generated/Segs_DFDX.f90');

%----------------------------------%
%     Parameter-Space Jacobian     %
%----------------------------------%
DFDP = jacobian(F_eqn, pvec);

% Write to file
fortran(DFDP, File='./fortran_generated/Segs_DFDP.f90');
