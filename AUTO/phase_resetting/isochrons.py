# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
import continuation_scripts.phase_reset.PTC_functions as phase_reset
# import plotting_scripts.phase_reset_plots as plot_PR
import plotting_scripts.isochron_plots as plot_iso

#=============================================================================#
#                  PHASE RESPONSE - PHASE TRANSITION CURVES                   #
#=============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

# %%
#-----------------------------------------------------------------------------#
#              Move Along Periodic Orbit (theta_old, theta_new)               #
#-----------------------------------------------------------------------------#
# We first move along with periodic orbit with zero-amplitude perturbation by
# freeing theta_old and theta_new. User defined points are saved at different
# values of theta_old for later isochron calculations. 

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run01_isochron_initial'

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

    # Read initial solution
    x1 = data_floquet['x1_read']
    x2 = data_floquet['x2_read']
    x3 = data_floquet['x3_read']

    #----------------------------#
    #     Initial Parameters     #
    #----------------------------#
    from numpy import pi

    # Integer for period
    k             = 50
    # \theta_old (where perturbation starts)
    theta_old     = 1.0
    # \theta_new (where segment comes back to \Gamma)
    theta_new     = 1.0
    # Distance from perturbed segment to \Gamma
    eta           = 0.0
    # Components of the perturbation vector
    d_x           = 0.0
    d_y           = 0.0
    d_z           = 0.0
    # Isochron components
    iso1          = x1[0].item()
    iso2          = x2[0].item()
    iso3          = x3[0].item()

    #----------------#
    #     Output     #
    #----------------#
    # Parameter vector
    p_out = { 1: gamma, 2: A, 3: B, 4: a,
              5: T, 6: k, 7: theta_old, 8: theta_new,
              9: mu_s, 10: eta,
             11: d_x, 12: d_y, 13: d_z,
             14: iso1, 15: iso2, 16: iso3}
    
    # Parameter names
    pnames_out = { 1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                   5: 'T', 6: 'k', 7: 'theta_old', 8: 'theta_new',
                   9: 'mu_s', 10: 'eta',
                  11: 'd_x', 12: 'd_y', 13: 'd_z',
                  14: 'iso1', 15: 'iso2', 16: 'iso3'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_PR, pnames_PR = set_parameters_PR()

# Calculate and write initial solution from previous run
phase_reset.write_initial_solution_phase_reset(par_PR[6])

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set saved points
from numpy import arange
SP_points = arange(0.0, 1.1, 0.1)
# SP_points = [1.0, 0.5]

# Copy continuation script
auto.copy('./continuation_scripts/phase_reset/', 'isochron_initial')

# Try set up phase reset calculation lol
run_new = auto.run(c='isochron_initial', PAR=par_PR, parnames=pnames_PR,
                   NTST=par_PR[6] * 50, JAC=1,
                   UZR={'theta_old': SP_points},
                   UZSTOP={'theta_old': 0.0})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_new, run_new_str)

#--------------#
#     Plot     #
#--------------#
label_plot = 1
# plot_PR.plot_phase_reset_phase_space(run_new, label_plot)
# plot_PR.plot_phase_reset_phase_space_2D(run_new, label_plot)

# Print clear line
print('\n')

# %%
#-----------------------------------------------------------------------------#
#                Second Continuation: Perturb Along x2 Direction              #
#-----------------------------------------------------------------------------#
# We pick an initial phase (theta_old) from the previous continuation, and 
# calculate the isochrons along the (x1, x2) plane, saving points along the
# x2 direction. We will pick up from these saved solutions in the next
# continuation to calculate isochrons in the (x1, x3) plane in slices along x2.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run02_isochron_dy_direction'
# Previous run name
run_old_str = 'run01_isochron_initial'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
# label_old = [sol['LAB'] for sol in run_old('UZ')]
# label_old = label_old[0]
label_old = 1

# Print to console
print('~~~ Phase Reset: Second Run ~~~')
print('Compute phase transition curve (PTC)')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set user defined points
SP_points = arange(-4.0, 0.03, 0.1)

# Continue from previous solution
run_new = auto.run(run_old(label_old), LAB=1,
                   DSMIN=1e-3, DS=1e-2, DSMAX=1e-1,
                   NMX=400, NPR=10,
                   ICP=['d_x', 'd_y', 'iso1', 'iso2', 'iso3', 'eta', 'mu_s', 'T'],
                   UZR={'d_y': SP_points},
                   UZSTOP={'d_x': [-5.0, 5.0], 'd_y': [-5.0, 5.0]})
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
# Plot single isochron
plot_iso.plot_single_isochron(run_new)

# # Plot phase space solution
# label_plot = 59
# plot_PR.plot_phase_reset_phase_space(run09_PTC_single, label_plot)

# Print clear line
print('\n')

# %%
#-----------------------------------------------------------------------------#
#       Third Continuation: Calculate Single Isochron in (x1, x3) Plane       #
#-----------------------------------------------------------------------------#
# We pick a point of the isochron along the x2-direction, and compute a slice 
# of the isochron in the (x1, x3) plane.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_single_isochron_plane'
# Previous run name
run_old_str = 'run02_isochron_dy_direction'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
# label_old = [sol['LAB'] for sol in run_old('UZ')]
# label_old = label_old[0]
label_old = 61

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
                   DSMIN=1e-2, DS=5e-1, DSMAX=1e0,
                   NMX=1000, NPR=50,
                   ICP=['d_x', 'd_z', 'iso1', 'iso2', 'iso3', 'eta', 'mu_s', 'T'])
# run_new += auto.run(DS='-')
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
# Plot single isochron
plot_iso.plot_single_isochron(run_new)

# # Plot phase space solution
# label_plot = 59
# plot_PR.plot_phase_reset_phase_space(run09_PTC_single, label_plot)

# Print clear line
print('\n')

# %%
#-----------------------------------------------------------------------------#
#     Fourth Continuation: Calculate Multiple Isochrons in (x1, x3) Plane     #
#-----------------------------------------------------------------------------#
# We compute along all saved points along the x2-direction, and compute a slice 
# of the isochron in the (x1, x3) plane.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run03_multi_isochron_plane'
# Previous run name
run_old_str = 'run02_isochron_dy_direction'
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
                        DSMIN=1e-2, DS=5e-1, DSMAX=1e0,
                        NMX=400, NPR=10,
                        ICP=['d_x', 'd_z', 'iso1', 'iso2', 'iso3', 'eta', 'mu_s', 'T'],
                        UZSTOP={'d_x': [-5.0, 5.0], 'd_z': [-5.0, 5.0]})
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
