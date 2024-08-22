%-------------------------------------------------------------------------%
%%                      Adjoint Equation Components                      %%
%-------------------------------------------------------------------------%
% State space variables
syms G Q I
xvec = [G; Q; I];
% Adjoint equation variables
syms w1 w2 w3
wvec = [w1; w2; w3];
% Parameters
syms gam A_pump B a
pvec = [gam; A_pump; B; a];

%-------------------------%
%     Field Equations     %
%-------------------------%
% Yamada model equations
F_vec = yamada(xvec, pvec);

% State space Jacobian
F_dfdx = jacobian(F_vec, xvec);

% Adjoint equations
A_vec = -transpose(F_dfdx) * wvec;

% Write to file
fortran(A_vec, File='./fortran_generated/adjoint.f90');

%------------------------------%
%     State-Space Jacobian     %
%------------------------------%
A_DFDX = jacobian(A_vec, [xvec; wvec]);

% Write to file
fortran(A_DFDX, File='./fortran_generated/adjoint_DFDX.f90');