# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
import continuation_scripts.phase_reset.PTC_functions as phase_reset
import plotting_scripts.PTC_plots as plot_PTC
import plotting_scripts.phase_reset_plots as plot_PR

#==============================================================================#
##                  PHASE RESPONSE - PHASE TRANSITION CURVES                  ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

# %%
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
run_new_str = 'run01_phase_reset_perturbation'

# Print to console
print('~~~ Phase Reset: First Run ~~~')
print('Continue in the perturbation ampltiude A_perturb')
print('Run name: {}'.format(run_new_str))

#-------------------#
#     Read Data     #
#-------------------#
# Define function for reading previous parameters and setting phase resetting
# parameters
def set_parameters_PR():
    from scipy.io import loadmat
    # Read Floquet data
    data_floquet = loadmat('./data_mat/floquet_solution.mat')
    
    # Read parameters
    gamma = data_floquet['gamma'].item()
    A     = data_floquet['A'].item()
    B     = data_floquet['B'].item()
    a     = data_floquet['a'].item()
    T     = data_floquet['T'].item()
    mu_s  = data_floquet['mu_s'].item()

    #----------------------------#
    #     Initial Parameters     #
    #----------------------------#
    from numpy import pi

    # Integer for period
    k             = 20
    # \theta_old (where perturbation starts)
    theta_old     = 1.0
    # \theta_new (where segment comes back to \Gamma)
    theta_new     = 1.0
    # Distance from perturbed segment to \Gamma
    eta           = 0.0
    # Size of perturbation
    A_perturb     = 0.0
    # Angle at which perturbation is applied?
    theta_perturb = 0.0
    # Azimuthal angle at which perturbation is applied?
    phi_perturb   = 0.0

    #----------------#
    #     Output     #
    #----------------#
    # Parameter vector
    p_out = { 1: gamma, 2: A, 3: B, 4: a,
              5: T, 6: k, 7: theta_old, 8: theta_new,
              9: mu_s, 10: eta,
             11: A_perturb, 12: theta_perturb, 13: phi_perturb}
    
    # Parameter names
    pnames_out = { 1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                   5: 'T', 6: 'k', 7: 'theta_old', 8: 'theta_new',
                   9: 'mu_s', 10: 'eta',
                  11: 'A_perturb', 12: 'theta_perturb', 13: 'phi_perturb'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_PR, pnames_PR = set_parameters_PR()

# Calculate and write initial solution from previous run
phase_reset.write_initial_solution_phase_reset(par_PR[6])

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import linspace, concatenate
SP_points = concatenate((linspace(0.0, 0.25, 25), linspace(0.30, 2.0, 25)))

# Copy continuation script
auto.copy('./continuation_scripts/phase_reset/', 'PTC_initial')

# Try set up phase reset calculation lol
run_new = auto.run(c='PTC_initial', PAR=par_PR, parnames=pnames_PR,
                   NMX=2000,
                   NTST=par_PR[6] * 50,
                   UZR={'A_perturb': SP_points},
                   UZSTOP={'A_perturb': max(SP_points) + 0.1})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

#--------------#
#     Plot     #
#--------------#
label_plot = 35
plot_PR.plot_phase_reset_phase_space(run_new, label_plot)
plot_PR.plot_phase_reset_phase_space_2D(run_new, label_plot)

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
run_new_str = 'run02_PTC_single'
# Previous run name
run_old_str = 'run01_phase_reset_perturbation'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
# label_old = [sol['LAB'] for sol in run_old('UZ')]
# label_old = label_old[0]
label_old = 32

# Print to console
print('~~~ Phase Reset: Second Run ~~~')
print('Compute phase transition curve (PTC)')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Continue from previous solution
run_new = auto.run(run_old(label_old), LAB=1,
                   ICP=['theta_old', 'theta_new', 'eta', 'mu_s', 'T'],
                   UZSTOP={'theta_old': [0.0, 2.0], 'theta_new': [-1.0, 12.0]},
                   # DSMIN=2.5e-2, DS=2.5e-2, DSMAX=2.5e-2,
                   NTST=100*35, NCOL=4, IAD=10,
                   ITMX=10, NWTN=3, ITNW=7,
                   NMX=8000, NPR=100,
                   JAC=0)
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
plot_PTC.plot_single_PTC(run_new)

# Plot phase space solution
label_plot = 59
plot_PR.plot_phase_reset_phase_space(run_new, label_plot)

# Plot temporal solution
plot_PR.plot_perturbed_states(run_new, label_plot)

# Plot perturbed intensity against theta_new
plot_PR.plot_theta_new_I_seg4(run_new)

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
run_new_str = 'run03_PTC_scan'
# Previous run name
run_old_str = 'run01_phase_reset_perturbation'
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
                        THL={'T': 0.0, 'mu_s': 0.0},
                        UZSTOP={'theta_old': [0.0, 2.0], 'theta_new': [-1.0, 12.0]},
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

#--------------#
#     Plot     #
#--------------#
# Save data
import python_files.save_PTC_scan_data as data_PTC
data_PTC.save_PTC_scan(run_new_str)

# Plot all PTC
# plot_PTC.plot_multi_PTC(run_new_str)


# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
