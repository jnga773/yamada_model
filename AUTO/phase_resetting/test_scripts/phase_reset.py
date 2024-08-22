# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import python_files.data_functions as data_funcs
import python_files.phase_reset_functions as phase_reset
import plotting_scripts.phase_reset_plots as plot_PR

#==============================================================================#
##                  PHASE RESPONSE - PHASE TRANSITION CURVES                  ##
#==============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

# %%
#------------------------------------------------------------------------------#
##              Move Along Periodic Orbit (theta_old, theta_new)              ##
#------------------------------------------------------------------------------#
# We first compute the phase response to a perturbation in a fixed direction
# applied at the zero-phase point, theta_old = 1.0 (or 0.0). We free  A_perturb
# and theta_old.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run08_PR_along_PO'

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
    from numpy import pi, cos, sin

    # Integer for period
    k             = 60
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

# Set saved points
from numpy import linspace
SP_points = linspace(1.0, 2.0, 21)

# Copy continuation script
auto.copy('./continuation_scripts/', 'phase_reset_initial')

# Try set up phase reset calculation lol
run08_PR_along_PO = auto.run(c='phase_reset_initial', PAR=par_PR, parnames=pnames_PR,
                             NTST=par_PR[6] * 50,
                             UZR={'theta_old': SP_points},
                             UZSTOP={'theta_old': max(SP_points) + 0.1})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run08_PR_along_PO, run_new_str)

#--------------#
#     Plot     #
#--------------#
# Plot label
label_plot = 11

# Plot phase space solution
plot_PR.plot_phase_reset_phase_space(run08_PR_along_PO, label_plot)
plot_PR.plot_phase_reset_phase_space_2D(run08_PR_along_PO, label_plot)
plot_PR.plot_perturbed_states(run08_PR_along_PO, label_plot)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
##                      Increase Perturbation Amplitude                       ##
#------------------------------------------------------------------------------#
# We fix the phase at which the perturbation is applied (theta_old) and increase
# the perturbation amplitude (A_perturb).

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run09_PR_increase_Ap'
# Previous run name
run_old_str = 'run08_PR_along_PO'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = 8

# Print to console
print('~~~ Phase Reset: Second Run ~~~')
print('Compute phase transition curve (PTC)')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from previous solution
run09_PR_increase_Ap = auto.run(run_old(label_old), LAB=1,
                            ICP=['A_perturb', 'theta_new', 'eta', 'mu_s', 'T'],
                            UZSTOP={'A_perturb': [0.0, 3.0]},
                            NMX=1000, NPR=10)

#-------------------#
#     Save Data     #
#-------------------#
data_funcs.save_move_data(run09_PR_increase_Ap, run_new_str)

# %%
#--------------#
#     Plot     #
#--------------#
# Plot label
label_plot = 25

# Plot phase space solution
plot_PR.plot_phase_reset_phase_space(run09_PR_increase_Ap, label_plot)
plot_PR.plot_phase_reset_phase_space_2D(run09_PR_increase_Ap, label_plot)
# plot_PR.plot_perturbed_states(run09_PR_increase_Ap, label_plot)

# Print clear line
print('\n')

# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
