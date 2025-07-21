#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 2 11:40:10 2025

@author: jnga773
"""

# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
import plotting_scripts.initial_PO_plots as plot_PO

#==============================================================================#
##                            INITIAL CONTINUATION                            ##
#==============================================================================#
# We calculate the initial periodic orbit of the Yamada model.

# Set parameters
gamma = 0.1
A     = 6.6
B     = 5.8
a     = 1.8
p0 = {1: gamma, 2: A, 3: B, 4: a}

# Initial solution is the 'off' state
x0 = [A, B, 0]

# Parameters to save the periodic orbit for
gamma_PO = 3.5e-2
A_PO     = 7.4

# %%
#------------------------------------------------------------------------------#
#                           Compute Equilibrium Point                          #
#------------------------------------------------------------------------------#
# We compute and continue the equilibrium point of the model.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run01_initial_EP'

# Print to console
print('~~~ Initial Periodic Orbit: First Run (c.initial_EP) ~~~')
print('Initial continuation from some point x0')
print('Run name: {}'.format(run_new_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# # Copy continuation script
auto.copy('./continuation_scripts/', 'initial_EP')

# Run the first continuation from the initial equilibrium point
run_new = auto.run(x0, PAR=p0, c='initial_EP')

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                         Continue from Branching Point                        #
#------------------------------------------------------------------------------#
# Continue from the branching point until we detect a Hopf bifurcation.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run02_branching_point'
# Previous run name
run_old_str = 'run01_initial_EP'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Second Run ~~~')
print('Continue bifurcations from the branching point')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Continue from branching point in run01_initial_EP
run_new = auto.run(run_old('BP1'), ISW=-1)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                            Continue Hopf To A_PO                             #
#------------------------------------------------------------------------------#
# We continue the Hopf bifurcation, varying the 'A' parameter until A = A_PO

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_hopf'
# Previous run name
run_old_str = 'run02_branching_point'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('HB1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Third Run ~~~')
print('Follow Hopf birfucation until z=-0.8')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Continue from Hop bifurcation
run_new = auto.run(run_old(label_old), ISW=2,
                   ICP=['A', 'gamma'],
                   UZSTOP={'gamma': [0.0, 0.4], 'A': [5.0, 20.0]},
                   DSMIN=5e-3, DS=5e-3, DSMAX=5e-3,
                   NMX=500, NPR=50,
                   UZR={'A': A_PO})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                            Hopf to Periodic Orbit                            #
#------------------------------------------------------------------------------#
# We compute a family of periodic orbits originating off the Hopf bifurcation.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run04_hopf_to_PO'
# Previous run name
run_old_str = 'run03_hopf'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('GH1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Fourth Run ~~~')
print('Continue periodic orbits from the Hopf bifurcation')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Follow periodic orbits
run_new = auto.run(run_old(label_old), ISW=-1, IPS=2, LAB=1,
                   ICP=['gamma'],
                   NTST=50, DSMIN=1e-1, DS=1e-1, DSMAX=1e-1,
                   NMX=500, NPR=50,
                   UZR={'gamma': gamma_PO})

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
run_new_str = 'run05_initial_PO'
# Previous run name
run_old_str = 'run04_hopf_to_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Fifth Run (c.initial_PO) ~~~')
print('Continue periodic orbit with shifted phase condition')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_PO, p_PO, pnames_PO = data_funcs.calc_initial_solution_PO(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_PO')

# Run continuation
run_new = auto.run(x_init_PO, PAR=p_PO, parnames=pnames_PO,
                   c='initial_PO', LAB=1)

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
run_new_str = 'run06_floquet_mu'
# Previous run name
run_old_str = 'run05_initial_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = 1

# Print to console
print('~~~ Floquet Bundle: First Run (c.initial_VAR) ~~~')
print('Calculate Floquet bundle (mu) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Calculate initial solution from previous run
x_init_VAR, p_VAR, pnames_VAR = data_funcs.calc_initial_solution_VAR(run_old(label_old))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_VAR')

# Run continuation
run_new = auto.run(x_init_VAR, PAR=p_VAR, parnames=pnames_VAR,
                   c='initial_VAR', NPR=50,
                   DSMIN=1e-4, DS=1e-4, DSMAX=1e-3)

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
run_new_str = 'run07_floquet_wnorm'
# Previous run name
run_old_str = 'run06_floquet_mu'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Floquet Bundle: Second Run (c.initial_VAR) ~~~')
print('Calculate Floquet bundle (w_norm) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

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
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
