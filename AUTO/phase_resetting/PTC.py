#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:24:15 2025

@author: jnga773
"""

# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
import plotting_scripts.PTC_plots as plot_PTC
import plotting_scripts.initial_PO_plots as plot_PO
import plotting_scripts.phase_reset_plots as plot_PR

# %%
#==============================================================================#
##                            INITIAL CONTINUATION                            ##
#==============================================================================#
# We calculate the initial periodic orbit of the Yamada model.

# Set parameters
gamma_PO = 3.5e-2
A_PO     = 7.4
B        = 5.8
a        = 1.8

# Parameter vector
p0 = {1: gamma_PO, 2: A_PO, 3: B, 4: a}

# Initial solution is the 'off' state
x0 = [10, 10, 10]

# State-space variable names
unames = {1: 'x1', 2: 'x2', 3: 'x3'}

# %%
#------------------------------------------------------------------------------#
#                    Confirm ODE45 Periodic Orbit Solution                     #
#------------------------------------------------------------------------------#
# Calculate the periodic orbit using MATLAB's ode45 function.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run01_initial_PO_ode45'

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Initial Periodic Orbit: First Run');
print('Find new periodic orbit');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Continuation parameters : {}'.format('gamma, T, A'));
print('=====================================================================');

#----------------------------#
#     Calculate Solution     #
#----------------------------#
# Solve using scipy's solve_ivp
x_init_solve_ivp = data_funcs.calc_initial_solution_solve_ivp(x0, p0)

# Parameter names
pnames_PO = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a', 11: 'T'}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points
UZR = {'A': A_PO, 'gamma': gamma_PO}
# Continuation parameters
pcont = ['gamma', 'T', 'A']
# Parameter range
prange = {'A': [5.0, 20.0], 'gamma': [0.0, 0.4]}

# Run continuation
run_new = auto.run(x_init_solve_ivp, e='./functions/yamada', IPS=2, IRS=0,
                   NPAR=len(pnames_PO), PAR=p0, parnames=pnames_PO, NDIM=len(unames), unames=unames,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   JAC=1, NBC=0, NINT=0,
                   NTST=50, DSMIN=1e-1, DS=1e-1, DSMAX=1e-1,
                   NMX=500, NPR=50)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

#--------------#
#     Plot     #
#--------------#
plot_PO.plot_hopf_to_PO_solution(run_new('UZ1'))

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#            Shifted Periodic Orbit with Zero-Phase Point Condition            #
#------------------------------------------------------------------------------#
# Continuing from a specified point in the previous continuation, we shift the
# origin of the periodic orbit to match the 'zero-phase point' condition, and
# continue the periodic orbit as a new BVP.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run02_initial_PO'
# Previous run name
run_old_str = 'run01_initial_PO_ode45'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Initial Periodic Orbit: Second Run');
print('Continue periodic orbit with shifted phase condition');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('A, gamma'));
print('=====================================================================');

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_PO, p_PO, pnames_PO = data_funcs.calc_initial_solution_PO(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points
UZR = {'A': A_PO, 'gamma': gamma_PO}
# Continuation parameters
pcont = ['gamma', 'T']
# Parameter range
prange = {'gamma': [0.0, 0.4]}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Run continuation
run_new = auto.run(x_init_PO, e='./functions/yamada_VAR', IPS=4, IRS=0, LAB=1,
                   NPAR=len(pnames_PO), PAR=p_PO, parnames=pnames_PO,
                   NDIM=len(unames), unames=unames,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   JAC=0, NBC=4, NINT=0,
                   NMX=10, NPR=1,
                   DSMIN=1e-3, DS=-1e-3, DSMAX=1e-2,
                   NCOL=4, IAD=1, NTST=50)

#-------------------#
#     Save Data     #
#-------------------#
# Save solution to MATLAB .mat file
data_funcs.save_data_PO(run_new(1), './data_mat/solution_PO.mat')

# Save data
data_funcs.save_move_data(run_new, run_new_str)

#--------------#
#     Plot     #
#--------------#
# Plot
plot_PO.plot_initial_PO_solution(run_new('EP1'))
# plot_PO.plot_initial_PO()

# Print clear line
print('\n')

# %%
#==============================================================================#
##                 COMPUTE FLOQUET BUNDLE AT ZERO PHASE POINT                 ##
#==============================================================================#
# Here we compute the stable Floquet bundle of the periodic orbit, as well
# as the perpendicular vector, w.

#------------------------------------------------------------------------------#
#               Continue the Floquet multiplier until mu_s = 1.0               #
#------------------------------------------------------------------------------#
# We now add the adjoint function and Floquet boundary conditions to
# compute the adjoint (left or right idk) eigenvectors and eigenvalues.
# This will give us the perpendicular vector to the tangent of the periodic
# orbit. However, this will only be for the eigenvector corresponding to
# the eigenvalue \mu = 1. Hence, here we continue in \mu (mu_f) until
# mu_f = 1.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_floquet_mu'
# Previous run name
run_old_str = 'run02_initial_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = 1

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Floquet Bundle: First Run');
print('Calculate stable Floquet bundle eigenvalue');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('mu_s, w_norm'));
print('=====================================================================');

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_VAR, p_VAR, unames_VAR, pnames_VAR = data_funcs.calc_initial_solution_VAR(run_old(label_old))

# Parameter names
pnames_PO = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a', 11: 'T'}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Saved points
UZR = {'mu_s': 1.0}
# Continuation parameters
pcont = ['mu_s', 'w_norm', 'gamma']
# Parameter range
prange = {'mu_s': [0.0, 1.1]}

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Run continuation
run_new = auto.run(x_init_VAR, e='./functions/yamada_VAR', IRS=0, IPS=4, ISW=1,
                   NPAR=len(p_VAR), PAR=p_VAR, parnames=pnames_VAR,
                   NDIM=len(unames_VAR), unames=unames_VAR,
                   ICP=pcont, UZSTOP=prange, UZR=UZR,
                   JAC=0, NBC=8, NINT=0,
                   NTST=50, NCOL=4, IAD=1,
                   DSMIN=5e-4, DS=1e-3, DSMAX=1e-2,
                   NMX=300, NPR=10)
#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                   Continue the w-vector until w_norm = 1.0                   #
#------------------------------------------------------------------------------#
# Having found the solution (branching point 'BP') corresponding to
# \mu = 1, we can continue in the norm of the vector w (w_norm), until the
# norm is equal to zero. Then we will have the correct perpendicular
# vector.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run04_floquet_wnorm'
# Previous run name
run_old_str = 'run03_floquet_mu'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Floquet Bundle: Second Run');
print('Grow norm of stable Floquet bundle vector');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('mu_s, w_norm'));
print('=====================================================================');

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Run continuation
run_new = auto.run(run_old(label_old), LAB=1,
                   ISW=-1, ILP=0,
                   DSMIN=5e-2, DS=1e-1, DSMAX=5e-1,
                   NMX=500, NPR=100,
                   UZSTOP={'w_norm': [0.0, 1.1]}, UZR={'w_norm': 1.0})

#-------------------#
#     Save Data     #
#-------------------#
# Save solution to MATLAB .mat file
data_funcs.save_data_VAR(run_new('UZ1'), './data_mat/solution_VAR.mat')

# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#==============================================================================#
##                  PHASE RESPONSE - PHASE TRANSITION CURVES                  ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

#------------------------------------------------------------------------------#
##                 First Continuation: Perturbation Amplitude                 ##
#------------------------------------------------------------------------------#
# We first compute the phase response to a perturbation in a fixed direction
# applied at the zero-phase point, theta_old = 1.0 (or 0.0). We free  A_perturb
# and theta_old.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run05_phase_reset_perturbation'
# Previous run name
run_old_str = 'run04_floquet_wnorm'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Phase Reset: First Run');
print('Continue in the perturbation ampltiude A_perturb');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('A_perturb, theta_new, eta, mu_s, T'));
print('=====================================================================');

#-------------------#
#     Read Data     #
#-------------------#
# Set initial phase resetting parameters
# Periodicity
k = 30

# Perturbation direction (in units of 2 \pi)
theta_perturb = 0.0
# theta_perturb = 0.25

# Calculate initial solution
x_init_PR, p_PR, unames_PR, pnames_PR = \
    data_funcs.calc_initial_solution_PR(run_old(label_old), k, theta_perturb,
                                        filename_out='./data_mat/solution_PR.dat',
                                        PTC_or_isochron='PTC')

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import linspace, concatenate, unique

# Saved points for large scan of G perturbation
SP_points = concatenate((linspace(0.0, 0.15, 20),
                         linspace(0.15, 1.0, 25),
                         linspace(1.0, 1.3, 25),
                         linspace(1.3, 2.0, 20)))
SP_points = concatenate((SP_points, [0.05, 0.1, 0.15, 0.5432, 1.0, 1.5, 2.0]))
SP_points = unique(SP_points)


# Set saved points
SP_points = {'A_perturb': SP_points}

# Set continuation parameters
pcont = ['A_perturb', 'theta_new', 'eta', 'mu_s', 'T']
# Set continuation stop points
prange = {'A_perturb': max(SP_points['A_perturb']) + 0.1}

run_new = auto.run(dat='./data_mat/solution_PR.dat',
                   e='./functions/yamada_PR', IRS=0, IPS=4, ISW=1, ILP=0, ISP=0,
                   NPAR=len(p_PR), PAR=p_PR, parnames=pnames_PR,
                   NDIM=len(unames_PR), unames=unames_PR,
                   JAC=1, NBC=22, NINT=0,
                   ICP=pcont, UZSTOP=prange, UZR=SP_points,
                   NTST=k*50, NCOL=4, IAD=10,
                   DSMIN=1e-2, DS=5e-2, DSMAX=1e0,
                   NMX=2000, NPR=100)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

#--------------#
#     Plot     #
#--------------#
# label_plot = 'UZ3'
# plot_PR.plot_phase_reset_phase_space(run_new(label_plot), label_plot)
# plot_PR.plot_phase_reset_phase_space_2D(run_new(label_plot), label_plot)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
##            Second Continuation: Phase Transition Curve (Single)            ##
#------------------------------------------------------------------------------#
# We then fix the perturbation amplitude A_perturb and free theta_old and
# theta_new to compute the phase transition curve (PTC).

# This is a test run of a SINGLE PTC.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run06_PTC_single'
# Previous run name
run_old_str = 'run05_phase_reset_perturbation'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
# label_old = run_old('UZ12')['LAB']
label_old = [sol['LAB'] for sol in run_old('UZ')]
label_old = label_old[2]
# label_old = 38

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Phase Reset: Second Run');
print('Compute phase transition curve (PTC) for single perturbation');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('theta_old, theta_new, eta, mu_s, T'));
print('=====================================================================');

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
from numpy import arange

# Set saved solutions for theta_old
SP_points = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
             1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
theta_old_stop = [0.0, 2.0]
theta_new_stop = [-1.0, 3.0]

# Saved points
UZR = {'theta_old': SP_points}
# Continuation parameters
pcont = ['theta_old', 'theta_new', 'eta', 'mu_s', 'T']
# Set continuation stop points
prange = {'theta_old': theta_old_stop, 'theta_new': theta_new_stop}

# Run continuation
run_scan = auto.run(run_old(label_old), LAB=1,
                    ICP=pcont, UZSTOP=prange, UZR=UZR,
                    DSMIN=1e-3, DS=1e-1, DSMAX=1e0,
                    NMX=4000, NPR=100)
run_scan += auto.run(DS='-')
# Merge results
run_new = auto.merge(run_new)
# Relabel results
run_new = auto.relabel(run_new)

#-------------------#
#     Save Data     #
#-------------------#
data_funcs.save_move_data(run_new, run_new_str)

#--------------#
#     Plot     #
#--------------#
# Plot single PTC solution
plot_PTC.plot_single_PTC(run_new)

# # Plot phase space solution
# label_plot = 17
# plot_PR.plot_phase_reset_phase_space(run_new, label_plot)

# # Plot temporal solution
# plot_PR.plot_perturbed_states(run_new, label_plot)

# # Plot perturbed intensity against theta_new
# plot_PR.plot_theta_new_I_seg4(run_new)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
##           Second Continuation: Phase Transition Curve (Multiple)           ##
#------------------------------------------------------------------------------#
# We then fix the perturbation amplitude A_perturb and free theta_old and
# theta_new to compute the phase transition curve (PTC).

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run07_PTC_scan'
# Previous run name
run_old_str = 'run05_phase_reset_perturbation'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ')
label_old = [sol['LAB'] for sol in label_old[:-1]]

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Phase Reset: Second Run');
print('Compute phase transition curve (PTC) for multiple perturbations');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Continuation parameters : {}'.format('theta_old, theta_new, eta, mu_s, T'));
print('=====================================================================');

#--------------------------------------#
#     Define Continuation Function     #
#--------------------------------------#
# Define function for parallelising
def calculate_PTC(i):
    """
    Run PTC continuation run for label 'i' in run_old.
    """
    from numpy import arange
    
    # This label
    this_label = label_old[i]

    # Run string identifier
    this_run = 'sol_' + str(i+1).zfill(3)

    # Print run information
    print('Continuing from point {} in run: {}'.format(this_label, run_old_str))
    
    # Set saved solutions for theta_old
    SP_points = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
    theta_old_stop = [0.0, 2.0]
    theta_new_stop = [-1.0, 3.0]
    
    # Saved points
    UZR = {'theta_old': SP_points}
    # Continuation parameters
    pcont = ['theta_old', 'theta_new', 'eta', 'mu_s', 'T']
    # Set continuation stop points
    prange = {'theta_old': theta_old_stop, 'theta_new': theta_new_stop}

    # Run continuation
    run_scan = auto.run(run_old(this_label), LAB=1,
                        ICP=pcont, UZSTOP=prange, UZR=UZR,
                        DSMIN=1e-3, DS=2e-3, DSMAX=5e-2,
                        EPSL=1e-7, EPSU=1e-7, EPSS=1e-4,
                        NMX=8000, NPR=100)
    run_scan += auto.run(DS='-')
    
    #-------------------#
    #     Save Data     #
    #-------------------#
    # Append runs and save data
    run_scan = auto.relabel(run_scan)
    # run_scan = auto.merge(run_scan)
    
    # Save data
    data_funcs.save_move_data(run_scan, '{}/{}'.format(run_new_str, this_run))

    # Print new line
    print('\n')

#--------------------------------------------#
#     Run Continuation: Regular For Loop     #
#--------------------------------------------#
# Regular for loop run
for i in range(len(label_old)):
    # Run continuation
    calculate_PTC(i)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
import continuation_scripts.save_PTC_scan_data as data_PTC
data_PTC.save_PTC_scan(run_new_str)

#--------------#
#     Plot     #
#--------------#
# Plot all PTC
# plot_PTC.plot_multi_PTC(run_new_str)


# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
