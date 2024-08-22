# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import python_files.data_functions as data_funcs
import python_files.initial_PO_functions as initial_PO
import plotting_scripts.initial_plots as plot_PO

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

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_EP')

# Run the first continuation from the initial equilibrium point
run01_initial_EP = auto.run(x0, PAR=p0, c='initial_EP')

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run01_initial_EP, run_new_str)

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

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from branching point in run01_initial_EP
run02_branching_point = auto.run(run01_initial_EP('BP1'), ISW=-1)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run02_branching_point, run_new_str)

# Print clear line
print('\n')

# %%
#------------------------------------------------------------------------------#
#                           Continue Hopf To z = -0.8                          #
#------------------------------------------------------------------------------#
# We continue the Hopf bifurcation, varying the 'z' parameter until z = -0.8

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

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from Hop bifurcation
run03_hopf = auto.run(run_old(label_old), ISW=2,
                      ICP=['A', 'gamma'],
                      UZSTOP={'gamma': [0.0, 0.4], 'A': [5.0, 20.0]},
                      DSMIN=5e-3, DS=5e-3, DSMAX=5e-3,
                      NMX=500, NPR=50,
                      UZR={'A': 6.715})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run03_hopf, run_new_str)

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
# label_old = run_old('GH1')
label_old = run_old('UZ1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Fourth Run ~~~')
print('Continue periodic orbits from the Hopf bifurcation')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Follow periodic orbits
print('Continuing periodic orbits ...')
run04_hopf_to_PO = auto.run(run_old(label_old), ISW=-1, IPS=2, LAB=1,
                            ICP=['gamma'],
                            NTST=50, DSMIN=1e-1, DS=1e-1, DSMAX=1e-1,
                            NMX=450, NPR=50,
                            UZR={'gamma': 7.52e-2})

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run04_hopf_to_PO, run_new_str)

# Write periodic orbit data to run06_initial_solution.dat
initial_PO.write_initial_solution_PO(run04_hopf_to_PO('UZ2'))

#--------------#
#     Plot     #
#--------------#
plot_PO.plot_hopf_to_PO_solution(run04_hopf_to_PO('UZ2'))

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
label_old = run_old('UZ2')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Periodic Orbit: Fifth Run (c.initial_PO) ~~~')
print('Continue periodic orbit with shifted phase condition')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
def read_parameters_PO(sol_in):
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['PAR(11)']


    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_PO, pnames_PO = read_parameters_PO(run_old(label_old))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_PO')

# Run continuation
run05_initial_PO = auto.run(c='initial_PO', PAR=par_PO, parnames=pnames_PO)

# Save data
data_funcs.save_move_data(run05_initial_PO, run_new_str)

# Write shifted periodic data and 'zero-ed' variational data as
# initial solution to the variational problem to ./data/run07_initial_solution.dat'
initial_PO.write_initial_solution_floquet(run05_initial_PO('EP1'))

#--------------#
#     Plot     #
#--------------#
# Save solution to MATLAB .mat file
initial_PO.save_PO_data_matlab(run05_initial_PO('EP1'))

# Plot
plot_PO.plot_initial_PO()

# Print clear line
print('\n')

#==============================================================================#
##                 COMPUTE FLOQUET BUNDLE AT ZERO PHASE POINT                 ##
#==============================================================================#
# Here we compute the stable Floquet bundle of the periodic orbit, as well
# as the perpendicular vector, w.

# %%
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
label_old = run_old('EP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Floquet Bundle: First Run (c.floquet_variational) ~~~')
print('Calculate Floquet bundle (mu) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#-------------------#
#     Read Data     #
#-------------------#
# Define function for reading previous parameters
def read_parameters_VAR(sol_in):
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']

    # Initial parameter values
    mu_s  = 0.8
    wnorm = 0.0

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T,
             6: mu_s, 7: wnorm}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T',
                  6: 'mu_s', 7: 'w_norm'}
    
    return p_out, pnames_out

# Read parameters from previous run
par_VAR, pnames_VAR = read_parameters_VAR(run_old(label_old))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'floquet_variational')

# Run continuation
run06_floquet_mu = auto.run(c='floquet_variational', PAR=par_VAR, parnames=pnames_VAR)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run06_floquet_mu, run_new_str)

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
print('~~~ Floquet Bundle: Second Run (c.floquet_variational) ~~~')
print('Calculate Floquet bundle (w_norm) ')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Try run this bitch
run07_floquet_wnorm = auto.run(run_old(label_old), LAB=1,
                               ISW=-1, ILP=0,
                               DSMIN=5e-2, DS=5e-2, DSMAX=5e-2,
                               NMX=200, NPR=10,
                               UZSTOP={'w_norm': [0.0, 1.1]}, UZR={'w_norm': 1.0})

#-------------------#
#     Save Data     #
#-------------------#
# Save solution to MATLAB .mat file
initial_PO.save_floquet_data_matlab(run07_floquet_wnorm('UZ1'))

# Save data
data_funcs.save_move_data(run07_floquet_wnorm, run_new_str)

# Print clear line
print('\n')

# %%
#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
