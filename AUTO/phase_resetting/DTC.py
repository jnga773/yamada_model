# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
# import plotting_scripts.PTC_plots as plot_PTC
# import plotting_scripts.phase_reset_plots as plot_PR
import plotting_scripts.DTC_plots as plot_DTC

# %%
#==============================================================================#
##                  PHASE RESPONSE - PHASE TRANSITION CURVES                  ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

#------------------------------------------------------------------------------#
##                     First Continuation: Move theta_old                     ##
#------------------------------------------------------------------------------#
# We first shift the starting point of the perturbation (theta_old) to 
# theta_old = 0.339413, which is where a pertubation at \theta_p = 45deg will
# intersect with the stable manifold of q.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run01_PR_move_theta_old'
# Previous run name
run_old_str = 'run07_floquet_wnorm'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

# Print to console
print('~~~ Phase Reset: First Run ~~~')
print('Move along periodic orbit')
print('Run name: {}'.format(run_new_str))

#-------------------#
#     Read Data     #
#-------------------#
# Set initial phase resetting parameters
# Periodicity
k = 50

# Perturbation direction (in units of 2 \pi)
# theta_perturb = 0.0
theta_perturb = 0.25

# Calculate initial solution
x_init_PR, p_PR, pnames_PR = \
    data_funcs.calc_initial_solution_PR(run_old(label_old), k, theta_perturb,
                                        filename_out='./data_mat/solution_PR.dat',
                                        PTC_or_isochron='PTC')

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import linspace, concatenate, unique

# Saved points for large scan of G perturbation
SP_points = 0.339413

# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_DTC')

# Try set up phase reset calculation lol
run_new = auto.run(dat='./data_mat/solution_PR.dat', PAR=p_PR, parnames=pnames_PR,
                   c='initial_DTC',
                   ICP=['theta_old', 'theta_new', 'eta', 'mu_s', 'T'],
                   NMX=300, NTST=k * 60,
                   UZR={'theta_old': SP_points},
                   UZSTOP={'theta_old': [0.0, 1.0]})

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
##                Second Continuation: Perturbation Amplitude                 ##
#------------------------------------------------------------------------------#
# We then fix this point and increase the perturbation amplitude (A_perturb).
# We save three-values to be used later to computer DTCs.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run02_DTC_A_perturb'
# Previous run name
run_old_str = 'run01_PR_move_theta_old'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
# label_old = run_old('UZ12')['LAB']
label_old = [sol['LAB'] for sol in run_old('UZ')]
label_old = label_old[0]
# label_old = 38

# Print to console
print('~~~ Phase Reset: Second Run ~~~')
print('Increase perturbation amplitude (A_perturb)')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
from numpy import arange

# Set saved solutions for A_perturb
SP_points = [0.1, 0.763750, 10.0]

# Continue from previous solution
run_new = auto.run(run_old(label_old), LAB=1,
                   ICP=['A_perturb', 'theta_new', 'eta', 'mu_s', 'T'],
                   DS='-',
                   UZSTOP={'A_perturb': [0.0, 10.1]},
                   UZR={'A_perturb': SP_points},
                   NMX=300, NPR=10)

#-------------------#
#     Save Data     #
#-------------------#
data_funcs.save_move_data(run_new, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
##                Third Continuation: Perturbation Direction                  ##
#------------------------------------------------------------------------------#
# We then fix the perturbation amplitude A_perturb and free theta_old and
# theta_new to compute the phase transition curve (PTC).

# This is a test run of a SINGLE PTC.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_DTC_test_single'
# Previous run name
run_old_str = 'run02_DTC_A_perturb'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
# label_old = run_old('UZ12')['LAB']
label_old = [sol['LAB'] for sol in run_old('UZ')]
label_old = label_old[1]
# label_old = 38

# Print to console
print('~~~ Phase Reset: Third Run ~~~')
print('Compute direction transition curve (DTC)')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
from numpy import arange

# Continue from previous solution
run_new = auto.run(run_old(label_old), LAB=1,
                   ICP=['theta_perturb', 'theta_new', 'eta', 'mu_s', 'T'],
                   DSMIN=1e-4, DS=1e-1, DSMAX=1e0,
                   UZSTOP={'theta_perturb': [-1, 1]},
                   NMX=4000, NPR=100)
run_new += auto.run(DS='-')
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
plot_DTC.plot_single_DTC(run_new)

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
run_new_str = 'run09_PTC_scan'
# Previous run name
run_old_str = 'run08_phase_reset_perturbation'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ')
label_old = [sol['LAB'] for sol in label_old[:-1]]

# Print to console
print('~~~ Phase Reset: Second Run ~~~')
print('Compute phase transition curve (PTC)')
print('Run name: {}'.format(run_new_str))
print('Continuing from (Saved Points) in run: {}'.format(run_old_str))

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
    SP_points = arange(0.1, 2.0, 0.1)
    theta_old_stop = [0.0, 2.0]
    theta_new_stop = [-1.0, 3.0]

    # Run continuation
    run_scan = auto.run(run_old(this_label), LAB=1,
                        ICP=['theta_old', 'theta_new', 'eta', 'mu_s', 'T'],
                        UZSTOP={'theta_old': theta_old_stop, 'theta_new': theta_new_stop},
                        UZR={'theta_old': SP_points},
                        DSMIN=1e-2, DS=1e-1, DSMAX=1e0,
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
