# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import python_files.write_data_functions as data_funcs
import plotting_scripts.PTC_plots as plot_PTC

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
run_new_str = 'run08_phase_reset_perturbation'
# Previous run name
run_old_str = 'run07_floquet_wnorm'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

# Print to console
print('~~~ Phase Reset: First Run ~~~')
print('Continue in the perturbation ampltiude A_perturb')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Define function for reading previous parameters and setting phase resetting
# parameters
def set_parameters_PR(sol_in):
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']
    mu_s  = sol_in['mu_s']

    #----------------------------#
    #     Initial Parameters     #
    #----------------------------#
    from numpy import pi, cos, sin

    # Integer for period
    k             = 35
    # Distance from perturbed segment to \Gamma
    eta           = 0.0
    # \theta_old (where perturbation starts)
    theta_old     = 1.0
    # \theta_new (where segment comes back to \Gamma)
    theta_new     = 1.0
    # Size of perturbation
    A_perturb     = 0.0
    # Angle at which perturbation is applied?
    theta_perturb = 0.5 * pi
    # Azimuthal angle at which perturbation is applied?
    phi_perturb = 0.0

    #----------------#
    #     Output     #
    #----------------#
    # Parameter vector
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T, 6: k, 7: mu_s, 8: eta,
             9: theta_old, 10: theta_new,
             11: A_perturb, 12: theta_perturb, 13: phi_perturb}
    # Parameter names
    pnames_out = { 1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                   5: 'T', 6: 'k', 7: 'mu_s', 8: 'eta',
                   9: 'theta_old', 10: 'theta_new',
                  11: 'A_perturb', 12: 'theta_perturb', 13: 'phi_perturb'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_PR, pnames_PR = set_parameters_PR(run_old(label_old))

# Calculate and write initial solution from previous run
data_funcs.write_initial_solution_phase_reset(run_old(label_old), k_in=par_PR[6])

# Set saved points
from numpy import linspace, concatenate
SP_points = concatenate((linspace(0.0, 0.25, 51), linspace(0.30, 25.0, 51)))

# Copy continuation script
auto.copy('./continuation_scripts/', 'PTC_initial')

# Try set up phase reset calculation lol
run08_phase_reset_perturbation = auto.run(c='PTC_initial', PAR=par_PR, parnames=pnames_PR, DS='-',
                                          NTST=par_PR[6] * 50,
                                          UZR={'A_perturb': SP_points},
                                          UZSTOP={'A_perturb': 25.1})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run08_phase_reset_perturbation, run_new_str)

#--------------#
#     Plot     #
#--------------#

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
##            Second Continuation: Phase Transition Curve (Single)            ##
#------------------------------------------------------------------------------#
# We then fix the perturbation amplitude A_perturb and free theta_old and
# theta_new to compute the phase transition curve (PTC).

# This is a test run of a SINGLE PTC.

# #------------------#
# #     Run Name     #
# #------------------#
# # This run name
# run_new_str = 'run09_PTC_single'
# # Previous run name
# run_old_str = 'run08_phase_reset_perturbation'
# run_old = data_funcs.bd_read(run_old_str)
# # Previous solution label
# label_old = [sol['LAB'] for sol in run_old('UZ')]
# label_old = label_old[0]

# # Print to console
# print('~~~ Phase Reset: Second Run ~~~')
# print('Compute phase transition curve (PTC)')
# print('Run name: {}'.format(run_new_str))
# print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

# #--------------------------#
# #     Run Continuation     #
# #--------------------------#
# # Continue from previous solution
# run09_PTC_single = auto.run(run_old(label_old), LAB=1,
#                             ICP=['theta_old', 'theta_new', 'eta', 'mu_s', 'T'],
#                             UZSTOP={'theta_old': [0.0, 2.0], 'theta_new': [-1.0, 12.0]},
#                             DSMIN=2.5e-2, DS=2.5e-2, DSMAX=2.5e-2,
#                             NTST=100*35, NCOL=4, IAD=10,
#                             ITMX=10, NWTN=3, ITNW=7,
#                             NMX=8000, NPR=100)
# run09_PTC_single += auto.run(DS='-')
# # Merge results
# run09_PTC_single = auto.merge(run09_PTC_single)
# # Relabel results
# run09_PTC_single = auto.relabel(run09_PTC_single)

# #-------------------#
# #     Save Data     #
# #-------------------#
# data_funcs.save_move_data(run09_PTC_single, run_new_str)

# #--------------#
# #     Plot     #
# #--------------#
# # Plot single PTC solution
# plot_PTC.plot_single_PTC(run09_PTC_single)

# # Print clear line
# print('\n')

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
    # This label
    this_label = label_old[i]

    # Run string identifier
    this_run = 'sol_' + str(i+1).zfill(3)

    # Print run information
    print('Continuing from point {} in run: {}'.format(this_label, run_old_str))

    # Run continuation
    run_scan = auto.run(run_old(this_label), LAB=1,
                        ICP=['theta_old', 'theta_new', 'eta', 'mu_s', 'T'],
                        UZSTOP={'theta_old': [0.0, 2.0], 'theta_new': [-1.0, 12.0]},
                        DSMIN=1.0e-2, DS=5.0e-2, DSMAX=1.0e-1,
                        NMX=8000, NPR=500)
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

#--------------#
#     Plot     #
#--------------#
# Save data
import python_files.save_PTC_scan_data as data_PTC
data_PTC.save_PTC_scan(run_new_str)

# Plot all PTC
# plot_PTC.plot_multi_PTC(run_new_str)

#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
