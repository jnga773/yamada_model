% State space dimension
xdim = 3;

%-------------------------------------------------------------------------%
%%                      Phase Resetting Components                       %%
%-------------------------------------------------------------------------%
%-------------------------------%
%     Initial-Point Vectors     %
%-------------------------------%
% Segment 1 - x(0)
x0_seg1 = sym('x0_seg1', [xdim, 1]);
assume(x0_seg1, 'real');

% Segment 1 - w(0)
w0_seg1 = sym('w0_seg1', [xdim, 1]);
assume(w0_seg1, 'real');

% Segment 2 - x(0)
x0_seg2 = sym('x0_seg2', [xdim, 1]);
assume(x0_seg2, 'real');

% Segment 2 - w(0)
w0_seg2 = sym('w0_seg2', [xdim, 1]);
assume(w0_seg2, 'real');

% Segment 3 - x(0)
x0_seg3 = sym('x0_seg3', [xdim, 1]);
assume(x0_seg3, 'real');

% Segment 4 - x(0)
x0_seg4 = sym('x0_seg4', [xdim, 1]);
assume(x0_seg4, 'real');

%-----------------------------%
%     Final-Point Vectors     %
%-----------------------------%
% Segment 1 - x(1)
x1_seg1 = sym('x1_seg1', [xdim, 1]);
assume(x1_seg1, 'real');

% Segment 1 - w(1)
w1_seg1 = sym('w1_seg1', [xdim, 1]);
assume(w1_seg1, 'real');

% Segment 2 - x(1)
x1_seg2 = sym('x1_seg2', [xdim, 1]);
assume(x1_seg2, 'real');

% Segment 2 - w(1)
w1_seg2 = sym('w1_seg2', [xdim, 1]);
assume(w1_seg2, 'real');

% Segment 3 - x(1)
x1_seg3 = sym('x1_seg3', [xdim, 1]);
assume(x1_seg3, 'real');

% Segment 4 - x(1)
x1_seg4 = sym('x1_seg4', [xdim, 1]);
assume(x1_seg4, 'real');

%--------------------%
%     Parameters     %
%--------------------%
% System parameters
syms gam A_pump B a
p_sys = [gam; A_pump; B; a];
assume(p_sys, 'real');

% Phase resetting parameters
syms T k
syms theta_old theta_new
syms mu_s eta
syms A_perturb theta_perturb phi_perturb
p_PR = [T; k; theta_old; theta_new;
        mu_s; eta;
        A_perturb; theta_perturb; phi_perturb];
assume(p_PR, 'real');

% Combined vector
pvec = [p_sys; p_PR];

%-------------------------------------------------------------------------%
%%                      Boundary Condition Encoding                      %%
%-------------------------------------------------------------------------%
%--------------------------%
%     Segments 1 and 2     %
%--------------------------%
% Vector field
F_vec = yamada(x0_seg1, p_sys);

% Boundary Conditions - Segments 1 and 2
bcs_seg12_1   = x0_seg1 - x1_seg2;
bcs_seg12_2   = x1_seg1 - x0_seg2;
bcs_seg12_3   = F_vec(1);

% Adjoint Boundary Conditions - Segments 1 and 2
a_bcs_seg12_1 = w1_seg1 - w0_seg2;
a_bcs_seg12_2 = (mu_s * w0_seg1) - w1_seg2;
a_bcs_seg12_3 = norm(w1_seg2) - 1;

%-------------------%
%     Segment 3     %
%-------------------%
bcs_seg3 = x1_seg3 - x0_seg1;

%-------------------%
%     Segment 4     %
%-------------------%
% Perturbation vector
d_vec = [cos(theta_perturb);
         0.0;
         sin(theta_perturb)];

% Boundary Conditions - Segment 4
bcs_seg4_1 = x0_seg4 - x0_seg3 - (A_perturb * d_vec);
bcs_seg4_2 = dot(x1_seg4 - x0_seg2, w0_seg2);
% bcs_seg4_3 = (norm(x1_seg4 - x0_seg2) ^ 2) - eta;
diff_vec = x1_seg4 - x0_seg2;
bcs_seg4_3 = (diff_vec(1) ^ 2) + (diff_vec(2) ^ 2) + (diff_vec(3) ^ 2) - eta;

%-----------------------%
%     Write to File     %
%-----------------------%
% Boundary condition vector
bcs = [bcs_seg12_1;
       bcs_seg12_2;
       bcs_seg12_3;
       a_bcs_seg12_1;
       a_bcs_seg12_2;
       a_bcs_seg12_3;
       bcs_seg3;
       bcs_seg4_1;
       bcs_seg4_2;
       bcs_seg4_3];

% Write bcs to Fortran file
fortran(bcs, File='./fortran_generated/bcs_Segs.f90');

%-------------------------------------------------------------------------%
%%                               JACOBIANS                               %%
%-------------------------------------------------------------------------%
% Initial vectors
U0_vec = [x0_seg1; w0_seg1;
          x0_seg2; w0_seg2;
          x0_seg3;
          x0_seg4];

% Final vectors
U1_vec = [x1_seg1; w1_seg1;
          x1_seg2; w1_seg2;
          x1_seg3;
          x1_seg4];

%-----------------------------%
%     State Vectors: Total    %
%-----------------------------%
% Combined state vector
uvec = [U0_vec; U1_vec];

% Calculate Jacobian
bcs_du = jacobian(bcs, uvec);

% Write to Fortran file
fortran(bcs_du, File='./fortran_generated/bcs_Segs_DFDU.f90');

%--------------------%
%     Parameters     %
%--------------------%
% Calculate total Jacobian
bcs_dp = jacobian(bcs, pvec);

% Write to Fortran file
fortran(bcs_dp, File='./fortran_generated/bcs_Segs_DFDP.f90');
