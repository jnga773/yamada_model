# # Append Python path to add AUTO functions
# import sys
# sys.path.append('/Users/jnga773/auto/07p/python')
# sys.path.append('/Users/jnga773/auto/07p/python/auto')

# Load extra functions
import auto
import python_files.write_data_functions as data_funcs
import plotting_scripts.bifurcation_plots as plot_bif

#==============================================================================#
#                              YAMADA LASER MODEL                              #
#==============================================================================#
# We compute a family of bifurcations from a driven-pulsed laser, described
# by the Yamada model of coupled differential equations:
#                    G' = \gamma (A - G - G I) ,
#                    Q' = \gamma (B - Q - a Q I) ,
#                    I' = (G - Q - 1) I ,
# where G is the gain, Q is the absorption, and I is the intensity of the
# laser. The system is dependent on four parameters: the pump current on
# the gain, A (or A); the relative absoprtion, B and a; and the decay
# time of the gain, \gamma.

# Set parameters
gamma = 0.1
A     = 6.6
B     = 5.8
a     = 1.8
p0 = {1: gamma, 2: A, 3: B, 4: a}

# Initial solution is the 'off' state
x0 = [A, B, 0]

#------------------------------------------------------------------------------#
#                           Compute Equilibrium Point                          #
#------------------------------------------------------------------------------#
# We compute and continue the equilibrium point of the model.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'initial_EP'

# Print to console
print('~~~ Initial Equilibrium Point: First Run (c.initial_EP) ~~~')
print('Initial continuation from some point x0')
print('Run name: {}'.format(run_new_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Copy continuation script
auto.copy('./continuation_scripts/', 'initial_EP')

# Run the first continuation from the initial equilibrium point
run_initial_EP = auto.run(x0, PAR=p0, c='initial_EP')

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_initial_EP, run_new_str)

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
#                         Continue from Branching Point                        #
#------------------------------------------------------------------------------#
# Continue from the branching point until we detect a Hopf bifurcation.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'branching_point'
# Previous run name
run_old_str = 'initial_EP'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Initial Equilibrium Point: Second Run ~~~')
print('Continue bifurcations from the branching point')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from branching point in run01_initial_EP
run_branching_point = auto.run(run_old(label_old), ISW=-1, LAB=1)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_branching_point, run_new_str)

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
#                           Continue Hopf To z = -0.8                          #
#------------------------------------------------------------------------------#
# We continue the Hopf bifurcation, varying the 'z' parameter until z = -0.8

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'hopf_bifurcations'
# Previous run name
run_old_str = 'branching_point'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('HB1')
label_old = label_old['LAB']

# Print to console
print('~~~ Hopf Bifrucations ~~~')
print('Follow Hopf birfucations')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from Hop bifurcation. Add a defined point (UZR) for double limit
# cycle continuation.
run_hopf = auto.run(run_old(label_old), ISW=2, LAB=1,
                    ICP=['A', 'gamma'],
                    UZSTOP={'gamma': [0.0, 0.4], 'A': [5.0, 15.0]},
                    DSMIN=5e-3, DS=5e-3, DSMAX=5e-3,
                    NMX=2000, NPR=100,
                    UZR={'A': 6.715})
# Run in opposite direction
run_hopf += auto.run(UZSTOP={'gamma': 0.0, 'A': 6.8},
                     DSMIN=5e-3, DS=-5e-3, DSMAX=5e-3)

#-------------------#
#     Save Data     #
#-------------------#
# Merge and relabel
run_hopf = auto.relabel(run_hopf)
run_hopf = auto.merge(run_hopf)

# Save data
data_funcs.save_move_data(run_hopf, run_new_str)

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
#                           Saddle Bifurcations: A_S                           #
#------------------------------------------------------------------------------#
# Continue the saddle-node point from run2 with a two parameter continuation to
# find the saddle-node line at A_S.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'saddle_AS'
# Previous run name
run_old_str = 'branching_point'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('LP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Saddle: A_S ~~~')
print('Follow line of saddle bifurcations')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from Hop bifurcation. Add a defined point (UZR) for double limit
# cycle continuation.
run_saddle_AS = auto.run(run_old(label_old), ISW=2,
                         ICP=['gamma', 'A'],
                         NMX=50, NPR=1,
                         UZR={'gamma': [1e-8, 0.4]})
run_saddle_AS += auto.run(DS='-', NMX=10)

#-------------------#
#     Save Data     #
#-------------------#
# Merge and relabel
run_saddle_AS = auto.relabel(run_saddle_AS)
run_saddle_AS = auto.merge(run_saddle_AS)

# Save data
data_funcs.save_move_data(run_saddle_AS, run_new_str)

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
#                       Transcritical Bifurcations: A_T                        #
#------------------------------------------------------------------------------#
# Continue the line of transcritical bifurcations at A = B+1.

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'transcritical_AT'
# Previous run name
run_old_str = 'initial_EP'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('BP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Transcritical: A_T ~~~')
print('Continue transcritical bifurcations for A_T line')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from Hop bifurcation. Add a defined point (UZR) for double limit
# cycle continuation.
run_saddle_AT = auto.run(run_old(label_old), ISW=2, ICP=['gamma', 'A'])
run_saddle_AT += auto.run(DS='-')

#-------------------#
#     Save Data     #
#-------------------#
# Merge and relabel
run_saddle_AT = auto.relabel(run_saddle_AT)
run_saddle_AT = auto.merge(run_saddle_AT)

# Save data
data_funcs.save_move_data(run_saddle_AT, run_new_str)

# Print clear line
print('\n')

#==============================================================================#
#                DOUBLE LIMIT CYCLE (SADDLE OF PERIODIC ORBITS)                #
#==============================================================================#
# Compute and follow the double limit cycles (saddle-nodes of periodic
# orbits). 

#------------------------------------------------------------------------------#
#                       Solve for Initial Periodic Orbit                       #
#------------------------------------------------------------------------------#
# We continue the Hopf bifurcation, varying the 'z' parameter until z = -0.8

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'double_limit_cycle_initial'
# Previous run name
run_old_str = 'hopf_bifurcations'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('UZ2')
label_old = label_old['LAB']

# Print to console
print('~~~ Double Limit Cycle: First Run ~~~')
print('Continue periodic orbits from Hopf until hit we hit a saddle')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from Hop bifurcation. Add a defined point (UZR) for double limit
# cycle continuation.
run_double_limit_initial = auto.run(run_old(label_old), LAB=1,
                                    ISW=-1, IPS=2,
                                    ICP=['gamma'],
                                    ISP=3,
                                    DSMIN=1e-3, DS=-1e-2, DSMAX=1e-1,
                                    NPR=10, NMX=200)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_double_limit_initial, run_new_str)

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
#                         Follow Saddle Periodic Orbits                        #
#------------------------------------------------------------------------------#
# We continue the Hopf bifurcation, varying the 'z' parameter until z = -0.8

#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'double_limit_cycle_saddle'
# Previous run name
run_old_str = 'double_limit_cycle_initial'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('LP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Double Limit Cycle: Second Run ~~~')
print('Continue saddle bifurcation of periodic orbits')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
# Continue from Hop bifurcation. Add a defined point (UZR) for double limit
# cycle continuation.
run_double_limit_saddle = auto.run(run_old(label_old), ISW=2, IPS=2, ICP=['A', 'gamma'])
# Rerun for some reason
run_double_limit_saddle = auto.run(run_double_limit_saddle('EP1'), LAB=1,
                                   ICP=['A', 'gamma'],
                                   NMX=5000, NPR=100,
                                   DSMIN=1e-2, DS=5e-2, DSMAX=1e-1)
run_double_limit_saddle += auto.run(DS='-')

#-------------------#
#     Save Data     #
#-------------------#
# Merge and relabel
run_double_limit_saddle = auto.relabel(run_double_limit_saddle)
run_double_limit_saddle = auto.merge(run_double_limit_saddle)

# Save data
data_funcs.save_move_data(run_double_limit_saddle, run_new_str)

# Print clear line
print('\n')

#==============================================================================#
#                  HOMOCLINIC (LARGE PERIOD PO APPROXIMATION)                  #
#==============================================================================#
# We compute a family of homoclinic orbits emanating from a Hopf bifurcation.
# We approximate the homoclinic orbit as a periodic orbit with a very large
# period.

#------------------------------------------------------------------------------#
#                      Compute Family of Periodic Orbits                       #
#------------------------------------------------------------------------------#
#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'homoclinic_PO_from_hopf'
# Previous run name
run_old_str = 'branching_point'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('HB1')
label_old = label_old['LAB']

# Print to console
print('~~~ Approximate Homoclinic: First Run ~~~')
print('Compute family of periodic orbits from the Hopf bifurcation')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
run_homoclinic_PO_from_hopf = auto.run(run_old(label_old), ISW=-1, IPS=2, LAB=1,
                                       ICP=['PERIOD', 'A'],
                                       NMX=5000, NPR=500,
                                       NTST=50,
                                       DSMIN=1e0, DS=1e0, DSMAX=1e0)

#-------------------#
#     Save Data     #
#-------------------#
# Save data
data_funcs.save_move_data(run_homoclinic_PO_from_hopf, run_new_str)

# Print clear line
print('\n')

#------------------------------------------------------------------------------#
#                         Large Period Periodic Orbit                          #
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
#                     Follow Family of "Homoclinic Orbits"                     #
#------------------------------------------------------------------------------#
#------------------#
#     Run Name     #
#------------------#
# This run name
run_new_str = 'homoclinic_approx_PO'
# Previous run name
run_old_str = 'homoclinic_PO_from_hopf'
run_old = data_funcs.bd_read(run_old_str)
# Previous solution label
label_old = run_old('EP1')
label_old = label_old['LAB']

# Print to console
print('~~~ Approximate Homoclinic: Second Run ~~~')
print('Follow family of homoclinic orbit')
print('Run name: {}'.format(run_new_str))
print('Continuing from point {} in run: {}'.format(label_old, run_old_str))

#--------------------------#
#     Run Continuation     #
#--------------------------#
run_homoclinic_approx = auto.run(run_old(label_old), ISW=1, IPS=2, LAB=1,
                                 ICP=['A', 'gamma'],
                                 UZSTOP={'A': [5.0, 6.81], 'gamma': [0.0, 0.4]},
                                 NMX=500, NPR=10,
                                 DSMIN=1e-3, DS=5e-3, DSMAX=1e-2)
run_homoclinic_approx += auto.run(DS='-')

#-------------------#
#     Save Data     #
#-------------------#
# Merge and relabel
run_homoclinic_approx = auto.relabel(run_homoclinic_approx)
run_homoclinic_approx = auto.merge(run_homoclinic_approx)

# Save data
data_funcs.save_move_data(run_homoclinic_approx, run_new_str)

# Print clear line
print('\n')


#==============================================================================#
#                          HOMOCLINIC (LIN'S METHOD)                           #
#==============================================================================#


#==============================================================================#
#                          SAVE AND PLOT BIFURCATION                           #
#==============================================================================#
# Save the data
run_names = {'H': 'hopf_bifurcations',
             'S': 'saddle_AS',
             'T': 'transcritical_AT',
             'D': 'double_limit_cycle_saddle',
             'L': 'homoclinic_approx_PO'}
plot_bif.save_data_mat(run_names)

# Plot diagrams
plot_bif.plot_bifurcation_diagram()
plot_bif.plot_bifurcation_diagram_zoomed()

#==============================================================================#
#                                 END OF FILE                                  #
#==============================================================================#
data_funcs.clean_directories()
