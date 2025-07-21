#==============================================================================#
#                FUNCTIONS USED IN THE PHASE RESET CALCULATIONS                #
#==============================================================================#
# This Python module contains the following functions:#
# - calc_stationary_points            : Calculates the three stationary
#                                       points of the Yamada model.
#
# - write_initial_solution_PO         : Reads the periodic orbit data from
#                                       the AUTO solution and writes it to
#                                       ./data/[XXX].dat.
#
# - calc_initial_solution_PO          : Reads the periodic orbit data from
#                                       the AUTO solution, rotates it, and sets
#                                       it as the initial solution.
#
# - save_data_PO_MATLAB               : Saves initial periodic orbit data to a 
#                                       MATLAB .mat data structure.
#
# - calc_initial_solution_VAR         : Reads the data from the previous AUTO
#                                       solution, and sets the initial
#                                       conditions.
#
# - save_data_VAR_MATLAB              : Saves the Floquet variational problem
#                                       as a MATLAB .mat file.

#------------------------------------------------------------------------------#
def calc_stationary_points(p_in):
    """
    Calculates the three stationary points of the Yamada model from analytic
    expression.

    Input
    ----------
    p_in : float, array
        Equation parameters.
    
    Output
    ----------
    x0_out : float, array
        The 'off' state
    xpos_out : float, array
        The higher-amplitude intensity state.
    xneg_out : float, array
        The lower-amplitude intensity state.  
    """
    from numpy import sqrt

    #---------------#
    #     Input     #
    #---------------#
    # Read the parameters
    gamma = p_in[0]
    A     = p_in[1]
    B     = p_in[2]
    a     = p_in[3]

    #--------------------------#
    #     Calculate Things     #
    #--------------------------#
    # Components for stationary points
    temp1 = -B - 1 - a + (a * A);
    temp1 = temp1 / (2 * a);

    temp2 = (B + 1 + a - (a * A)) ** 2;
    temp2 = temp2 - 4 * (1 + B - A) * a;
    temp2 = sqrt(temp2);
    temp2 = temp2 / (2 * a);

    # The two non-trivial equilibrium values of intensity
    I_pos = temp1 + temp2;
    I_neg = temp1 - temp2;

    #----------------#
    #     Output     #
    #----------------#
    x0_out   = [A, B, 0]
    xpos_out = [A / (1 + I_pos), B / (1 + (a * I_pos)), I_pos]
    xneg_out = [A / (1 + I_neg), B / (1 + (a * I_neg)), I_neg]

    return x0_out, xpos_out, xneg_out

#------------------------------------------------------------------------------#
def calc_initial_solution_PO(sol_in, filename_out):
    """
    Reads the periodic orbit data from [sol_in], rotates the periodic orbit,
    and writes the initial time and state space solution to [filename_out]

    Input
    ----------
    sol_in : AUTO-continuation solution
        AUTO generated solution for the periodic orbit (like run('UZ')).
        Probably run05_initial_PO('UZ1'), along with the stationary points.
    filename_out : str
        Filename of the .dat file to save the data to.
    
    Output
    ----------
    x_init_out: array
        Initial solution as columns [time, state-space solution].
    p_out: array
        Array of the initial parameter values
    pnames_out: cell
        Cell containing the parameter names and corresponding indices.
    """
    from numpy import argmax, concatenate, array
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t_read = sol_in['t']

    # Read data
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['PAR(11)']

    # Calculate stationary points
    x0, xpos, xneg = calc_stationary_points([gamma, A, B, a])

    #--------------------#
    #     Shift Data     #
    #--------------------#
    # Find the point where the first component is the maximum
    max_idx = argmax(x1_read)

    # Shift data
    if max_idx > 0:
        # Shift time array
        t = [t_read[max_idx:] - t_read[max_idx],
             t_read[1:max_idx+1] + (t_read[-1] - t_read[max_idx])]
        
        # Shift around arrays
        x1 = [x1_read[max_idx:], x1_read[1:max_idx+1]]
        x2 = [x2_read[max_idx:], x2_read[1:max_idx+1]]
        x3 = [x3_read[max_idx:], x3_read[1:max_idx+1]]

        # Concatenate arrays
        t = concatenate((t[0], t[1]))
        x1 = concatenate((x1[0], x1[1]))
        x2 = concatenate((x2[0], x2[1]))
        x3 = concatenate((x3[0], x3[1]))

    #-----------------------#
    #     Write to file     #
    #-----------------------#
    with open(filename_out, 'w') as file_data:
        for i in range(len(t)):
            str_write = ("{:>20.15f} \t "
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t[i],
                                  x1[i], x2[i], x3[i],
                                  x0[0], x0[1], x0[2],
                                  xneg[0], xneg[1], xneg[2],
                                  xpos[0], xpos[1], xpos[2])
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

    #----------------#
    #     Output     #
    #----------------#
    # State space solution
    x_init_out = [t, x1, x2, x3]
    x_init_out = array(x_init_out)

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T}
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T'}

    return x_init_out, p_out, pnames_out

#------------------------------------------------------------------------------#
def save_data_PO_MATLAB(sol_in, filename_out):
    """
    Writes the periodic orbit and stationary point solution to a MATLAB .mat
    file.

    Input
    -------
    sol_in : AUTO-continuation solution
        AUTO generated solution for the periodic orbit (like run('UZ')).
        Probably run05_initial_PO('UZ1'), along with the stationary points.
    filename_out : str
        Filename of the .mat file to save the data to.
    """
    from scipy.io import savemat
    from numpy import array

    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    #------------------------#
    #     Periodic Orbit     #
    #------------------------#
    # Get state space data
    x1 = sol_in['x1']
    x2 = sol_in['x2']
    x3 = sol_in['x3']

    # State vector
    xbp_PO = array([x1, x2, x3]).T
    
    #----------------------------#
    #     Equilibrium Points     #
    #----------------------------#
    # x0
    ss_1 = sol_in['x0_1']
    ss_2 = sol_in['x0_2']
    ss_3 = sol_in['x0_3']
    # State vector
    x0 = array([ss_1, ss_2, ss_3])

    # xpos
    ss_1 = sol_in['xpos_1']
    ss_2 = sol_in['xpos_2']
    ss_3 = sol_in['xpos_3']
    # State vector
    xpos = array([ss_1, ss_2, ss_3])

    # xneg
    ss_1 = sol_in['xneg_1']
    ss_2 = sol_in['xneg_2']
    ss_3 = sol_in['xneg_3']
    # State vector
    xneg = array([ss_1, ss_2, ss_3])

    #--------------------#
    #     Parameters     #
    #--------------------#
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']

    # Time data
    tbp = sol_in['t']
    # Multiply by period
    tbp *= T

    #-------------------#
    #     Save Data     #
    #-------------------#
    # Dictionary for data
    mat_out = {'tbp': tbp, 'xbp_PO': xbp_PO,
               'x0': x0, 'xneg': xneg, 'xpos': xpos,
               'p0_PO': [gamma, A, B, a]}
    
    # Save data
    savemat(filename_out, mat_out) 

#------------------------------------------------------------------------------#
def calc_initial_solution_VAR(sol_in, filename_out):
    """
    Reads the periodic orbit data from [sol_in], rotates the periodic orbit,
    and writes the initial time and state space solution to [filename_out]

    Input
    ----------
    sol_in : AUTO-continuation solution
        AUTO generated solution for the rotated periodic orbit and variational 
        problem solution. Probably run07_floquet_wnorm('UZ1').
    filename_out : str
        Filename of the .dat file to save the data to.
    
    Output
    ----------
    x_init_out: array
        Initial solution as columns [time, state-space solution].
    p_out: array
        Array of the initial parameter values
    pnames_out: cell
        Cell containing the parameter names and corresponding indices.
    """
    from numpy import zeros, array

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t = sol_in['t']

    # Read data
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    # Zeros for perpindicular solution
    w1 = zeros(shape=(len(x1_read)))
    w2 = zeros(shape=(len(x1_read)))
    w3 = zeros(shape=(len(x1_read)))

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']

    # Initial parameter values
    mu_s  = 0.8
    wnorm = 0.0

    #-----------------------#
    #     Write to file     #
    #-----------------------#
    with open(filename_out, 'w') as file_data:
        for i in range(len(t)):
            str_write = ("{:>20.15f} \t "
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t[i],
                                  x1_read[i], x2_read[i], x3_read[i],
                                  0.0, 0.0, 0.0)
            
            # Write to file
            file_data.write(str_write)

    #----------------#
    #     Output     #
    #----------------#
    # State space solution
    x_init_out = [t, x1_read, x2_read, x3_read, w1, w2, w3]
    x_init_out = array(x_init_out, dtype='float')

    # Parameter values
    p_out = {1: gamma, 2: A, 3: B, 4: a,
             5: T,
             6: mu_s, 7: wnorm}
    
    # Parameter names
    pnames_out = {1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                  5: 'T',
                  6: 'mu_s', 7: 'w_norm'}

    return x_init_out, p_out, pnames_out

#------------------------------------------------------------------------------#
def save_VAR_data_matlab(bd_data_in):
    """
    Read the final solution from run07_floquet_wnorm.
    """
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Grab the periodic orbit data file
    sol = bd_data_in

    # Read time data
    t_read = sol['t']

    # Periodic orbit solution
    x1_read = sol['x1']
    x2_read = sol['x2']
    x3_read = sol['x3']

    # Perpindicular vector solution
    wn_1_read = sol['wn_1']
    wn_2_read = sol['wn_2']
    wn_3_read = sol['wn_3']

    #-----------------------#
    #     Write to file     #
    #-----------------------#
    filename_data = './data/VAR_data_test.dat'
    with open(filename_data, 'w') as file_data:
        for i in range(len(t_read)):
            str_write = ("{:>20.15f} \t "
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t"
                         "{:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
                         ).format(t_read[i],
                                  x1_read[i], x2_read[i], x3_read[i],
                                  wn_1_read[i], wn_2_read[i], wn_3_read[i])
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

#------------------------------------------------------------------------------#
def save_data_VAR_MATLAB(sol_in, filename_out):
    """
    Writes the periodic orbit and Floquet variation problem data to a MATLAB
    .mat file/

    Input
    -------
    sol_in : AUTO-continuation solution
        AUTO generated solution for the periodic orbit (like run('UZ')).
        Probably run05_initial_PO('UZ1'), along with the stationary points.
    filename_out : str
        Filename of the .mat file to save the data to.
    """
    from scipy.io import savemat
    from numpy import array

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t_read = sol_in['t']

    #--------------------#
    #     State data     #
    #--------------------#
    # Periodic orbit solution
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    # Variational vector solution
    wn1_read = sol_in['wn_1']
    wn2_read = sol_in['wn_2']
    wn3_read = sol_in['wn_3']

    #--------------------#
    #     Parameters     #
    #--------------------#
    # Read parameters
    gamma = sol_in['gamma']
    A     = sol_in['A']
    B     = sol_in['B']
    a     = sol_in['a']
    T     = sol_in['T']
    mu_s  = sol_in['mu_s']

    #-------------------#
    #     Save Data     #
    #-------------------#
    # Dictionary for data
    mat_out = {'t_read': t_read,
               'x1_read': x1_read, 'x2_read': x2_read, 'x3_read': x3_read,
               'wn1_read': wn1_read, 'wn2_read': wn2_read, 'wn3_read': wn3_read,
               'gamma': gamma, 'A': A, 'B': B, 'a': a, 'T': T, 'mu_s': mu_s}
    
    # Save data
    savemat(filename_out, mat_out, oned_as='column') 
