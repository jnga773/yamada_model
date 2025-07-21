# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import continuation_scripts.data_functions as data_funcs
import continuation_scripts.phase_reset as phase_reset
import plotting_scripts.phase_reset_plots as plot_PR
import plotting_scripts.isochron_plots as plot_iso

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
#=============================================================================#
#                         PHASE RESPONSE - ISOCHRONES                         #
#=============================================================================#
# We set up the phase resetting problem by creating four segments of the
# periodic orbit, with boundary conditions described in a paper somewhere.

#------------------------------------------------------------------------------#
##               First Continuation: Follow Along Periodic Orbit              ##
#------------------------------------------------------------------------------#
# We first compute the phase response to a perturbation of amplitude zero along
# the periodic orbit, saving at different theta_old phase points.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run08_isochron_initial'
# Previous run name
run_old_str = 'run07_floquet_wnorm'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ1')
label_old = label_old['LAB']

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Isochron: First Run');
print('Change phase along periodic orbit');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('theta_old, theta_new, eta, mu_s, T'));
print('=====================================================================');

#-------------------#
#     Read Data     #
#-------------------#
# Set initial phase resetting parameters
# Periodicity
k = 35

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
from numpy import arange
# SP_points = arange(0.0, 1.1, 0.1)
SP_points = [0.1858, 0.6768, 1.0]
# SP_points = [1.0, 0.5]

# Copy continuation script
auto.copy('./continuation_scripts/', 'isochron_initial')

# Try set up phase reset calculation lol
run_new = auto.run(dat='./initial_solutions/PR.dat', PAR=p_PR, parnames=pnames_PR,
                   c='isochron_initial',
                   NTST=k * 50,
                   UZR={'theta_old': SP_points},
                   UZSTOP={'theta_old': 0.0})

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
run_new_str = 'run09_isochron_dy_direction'
# Previous run name
run_old_str = 'run08_isochron_initial'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = [sol['LAB'] for sol in run_old('UZ')]
label_old = label_old[1]
label_old = 1

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Isochron: Second Run');
print('Fix theta_old and theta_new, and follow isochrons in one direction');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('d_x, d_y, eta, mu_s, iso1, iso2, iso3'));
print('=====================================================================');

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Set user defined points
SP_points = arange(0.0, 4.2, 0.2)

# Continue from previous solution
run_new = auto.run(run_old(label_old), LAB=1,
                   DSMIN=1e-3, DS=1e-2, DSMAX=1e-1,
                   NMX=600, NPR=50,
                   ICP=['d_x', 'd_y', 'iso1', 'iso2', 'iso3', 'eta', 'mu_s', 'T'],
                   UZR={'iso2': SP_points},
                   UZSTOP={'iso1': [-4.0, 5.0], 'iso2': [-2.0, 5.0]})

# run_new += auto.run(DS='-')
# # Merge results
# run_new = auto.merge(run_new)
# # Relabel results
# run_new = auto.relabel(run_new)

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
run_new_str = 'run10_single_isochron_plane'
# Previous run name
run_old_str = 'run09_isochron_dy_direction'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = [sol['LAB'] for sol in run_old('UZ')]
label_old = label_old[5]

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Isochron: Third Run');
print('Compute isochron slice in another plane');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('d_x, d_z, eta, mu_s, iso1, iso2, iso3'));
print('=====================================================================');

#-------------------------------#
#     Run AUTO Continuation     #
#-------------------------------#
# Continue from previous solution
run_new = auto.run(run_old(label_old), LAB=1,
                   NMX=1500, NPR=250,
                   ICP=['d_x', 'd_z', 'iso1', 'iso2', 'iso3', 'eta', 'mu_s', 'T'],
                   UZSTOP={'iso1': [-4.0, 6.0], 'iso2': [-4.0, 6.0]})
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
#     Fourth Continuation: Calculate Multiple Isochrons in (x1, x3) Plane     #
#-----------------------------------------------------------------------------#
# We compute along all saved points along the x2-direction, and compute a slice 
# of the isochron in the (x1, x3) plane.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'run10_multi_isochron_plane'
# Previous run name
run_old_str = 'run09_isochron_dy_direction'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = [sol['LAB'] for sol in run_old('UZ')]

#--------------------------#
#     Print to Console     #
#--------------------------#
print('=====================================================================');
print('Isochron: Fourth Run');
print('Calculate slices of the isochrons');
print('---------------------------------------------------------------------');
print('This run name           : {}'.format(run_new_str));
print('Previous run name       : {}'.format(run_old_str));
print('Previous label_solution : {}'.format(label_old));
print('Continuation parameters : {}'.format('d_x, d_z, eta, mu_s, iso1, iso2, iso3'));
print('=====================================================================');

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
                        # DSMIN=1e-2, DS=5e-1, DSMAX=1e0,
                        NMX=2000, NPR=250,
                        ICP=['d_x', 'd_z', 'iso1', 'iso2', 'iso3', 'eta', 'mu_s', 'T'],
                        UZSTOP={'iso1': [-4.0, 6.0], 'iso2': [-4.0, 6.0]})
    run_scan += auto.run(DS='-')
    
    #-------------------#
    #     Save Data     #
    #-------------------#
    # Merge run
    run_scan = auto.merge(run_scan)
    # Append runs and save data
    run_scan = auto.relabel(run_scan)
    
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
import continuation_scripts.phase_reset.save_isochron_data as save_iso
save_iso.save_isochron_scan_data(run_new_str)

# Plot all isochrons
# plot_iso.plot_single_isochron(run_new)
plot_iso.plot_multi_isochron(run_new_str)


# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
