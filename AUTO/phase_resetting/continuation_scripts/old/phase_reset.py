#==============================================================================#
#                FUNCTIONS USED IN THE PHASE RESET CALCULATIONS                #
#==============================================================================#
# This Python module contains the following functions:
#
# - calc_initial_solution_PR          : Calculates and sets the initial
#                                       solution to solve for the phase
#                                       resetting problems, and write it
#                                       to ./data/[XXX].dat file.

#------------------------------------------------------------------------------#
# Calculate the initial solutions to the phase resetting problem.
def calc_initial_solution_PR(sol_in, k_in, theta_in, PTC_or_ISO, filename_out):
    """
    Reads the variational problem solution from [sol_in], computes the initial
    solution to the phase resetting problem, and writes it to a .dat file
    [filename_out]. Input the periodicity [k_in] and the direction of the
    perturbation in the G-I plane [theta_in], with vector:
            d = (cos(theta_in) , 0 , sin(theta_in)).

    Input
    ----------
    sol_in : AUTO-continuation solution
        AUTO generated solution for the periodic orbit (like run('UZ')).
        Probably run05_initial_PO('UZ1'), along with the stationary points.
    k_in : int
        Periodicity of the solution.
    theta_in : float
        Perturbation directional angle in the G-I plane.
    PTC_or_ISO : str {'PTC', 'ISO'}
        The input parameter for the PTC and isochron calculations are different
        so this will determine which one we want.
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
    from numpy import concatenate, interp, ones, array

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read time data
    t_read = sol_in['t']

    # Periodic orbit solution
    x1_read = sol_in['x1']
    x2_read = sol_in['x2']
    x3_read = sol_in['x3']

    # Perpindicular vector solution
    wn1_read = sol_in['wn_1']
    wn2_read = sol_in['wn_2']
    wn3_read = sol_in['wn_3']

    # Initial points
    x1_0 = x1_read[0]
    x2_0 = x2_read[0]
    x3_0 = x3_read[0]
    wn1_0 = wn1_read[0]
    wn2_0 = wn2_read[0]
    wn3_0 = wn3_read[0]

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
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
    # Integer for period
    k             = k_in
    # \theta_old (where perturbation starts)
    theta_old     = 1.0
    # \theta_new (where segment comes back to \Gamma)
    theta_new     = 1.0
    # Distance from perturbed segment to \Gamma
    eta           = 0.0
    # Size of perturbation
    A_perturb     = 0.0
    # Angle at which perturbation is applied?
    theta_perturb = theta_in
    # Azimuthal angle at which perturbation is applied?
    phi_perturb   = 0.0
    # For the isochron calculations, perturbation vector components
    d_x           = 0.0
    d_y           = 0.0
    d_z           = 0.0
    
    if PTC_or_ISO == 'ISO':
        # Isochron components
        iso1          = x1_0
        iso2          = x2_0
        iso3          = x3_0

    #--------------------------#
    #     "Periodise" Data     #
    #--------------------------#
    # Set initial time array
    t_period  = t_read

    # Set inintial state date
    x1_period = x1_read
    x2_period = x2_read
    x3_period = x3_read

    # If k > 1, cycle through and keep appending solutions
    if k_in > 1:
        # Cycle through number of periods
        for i in range(k_in-1):
            # Append time
            t_period = concatenate((t_read, 1 + t_period[1:]))
            
            # Segment 4: State space data
            x1_period = concatenate((x1_read, x1_period[1:]))
            x2_period = concatenate((x2_read, x2_period[1:]))
            x3_period = concatenate((x3_read, x3_period[1:]))
            
        # Normalise time data
        t_period *= 1 / k_in

    #--------------------------#
    #     Interpolate Data     #
    #--------------------------#
    # Interpolate: State space
    x1_interp  = interp(t_period, t_read, x1_read)
    x2_interp  = interp(t_period, t_read, x2_read)
    x3_interp  = interp(t_period, t_read, x3_read)

    # Interpolate: Perpindicular space
    wn1_interp = interp(t_period, t_read, wn1_read)
    wn2_interp = interp(t_period, t_read, wn2_read)
    wn3_interp = interp(t_period, t_read, wn3_read)

    # New ones matrix
    ones_mat = ones(shape=(len(t_period)))

    # Multiply through single values
    x1_extended  = x1_0 * ones_mat
    x2_extended  = x2_0 * ones_mat
    x3_extended  = x3_0 * ones_mat

    wn1_extended = wn1_0 * ones_mat
    wn2_extended = wn2_0 * ones_mat
    wn3_extended = wn3_0 * ones_mat
    
    #--------------------#
    #     Write Data     #
    #--------------------#
    # Time
    t_seg = t_period

    # Segment 1
    x_seg1 = [x1_interp, x2_interp, x3_interp,
              wn1_interp, wn2_interp, wn3_interp]
    x_seg1 = array(x_seg1).T

    # Segment 2
    x_seg2 = [x1_extended, x2_extended, x3_extended,
              wn1_extended, wn2_extended, wn3_extended]
    x_seg2 = array(x_seg2).T

    # Segment 3
    x_seg3 = [x1_extended, x2_extended, x3_extended]
    x_seg3 = array(x_seg3).T

    # Segment 4
    x_seg4 = [x1_period, x2_period, x3_period]
    x_seg4 = array(x_seg4).T

    with open(filename_out, 'w') as file_data:
        for i in range(len(t_seg)):
            # Break up each column section into separate pieces

            # Time
            t_write    = "{:>15.10f}".format(t_seg[i])
            # Segment 1
            seg1_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f} \t"
                          "{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg1[i, 0], x_seg1[i, 1], x_seg1[i, 2],
                                   x_seg1[i, 3], x_seg1[i, 4], x_seg1[i, 5])
            # Segment 2
            seg2_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f} \t"
                          "{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg2[i, 0], x_seg2[i, 1], x_seg2[i, 2],
                                   x_seg2[i, 3], x_seg2[i, 4], x_seg2[i, 5])
            # Segment 3
            seg3_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg3[i, 0], x_seg3[i, 1], x_seg3[i, 2])
            # Segment 4
            seg4_write = ("{:>15.10f} \t {:>15.10f} \t {:>15.10f}"
                          ).format(x_seg4[i, 0], x_seg4[i, 1], x_seg4[i, 2])

            # Total string
            str_write = ("{} \t {} \t {} \t {} \t {} \n"
                         ).format(t_write, seg1_write, seg2_write, seg3_write, seg4_write)
            
            # Write to file
            file_data.write(str_write)

    # Close file
    file_data.close()

    #----------------#
    #     Output     #
    #----------------#
    # Initial vector
    x_init_out = [t_period,
                  # Segment '1'
                  x1_interp, x2_interp, x3_interp,
                  wn1_interp, wn2_interp, wn3_interp,
                  # Segment '2'
                  x1_extended, x2_extended, x3_extended, 
                  wn1_extended, wn2_extended, wn3_extended,
                  # Segment '3'
                  x1_extended, x2_extended, x3_extended,
                  # Segment '4'
                  x1_period, x2_period, x3_period]
    x_init_out = array(x_init_out, dtype='float')

    # Parameter vector
    p_out = { 1: gamma, 2: A, 3: B, 4: a,
            5: T, 6: k, 7: theta_old, 8: theta_new,
            9: mu_s, 10: eta}
    
    # Parameter names
    pnames_out = { 1: 'gamma', 2: 'A', 3: 'B', 4: 'a',
                   5: 'T', 6: 'k', 7: 'theta_old', 8: 'theta_new',
                   9: 'mu_s', 10: 'eta'}

    if PTC_or_ISO == 'PTC':
        # Parameter vector
        p_out[11] = A_perturb
        p_out[12] = theta_perturb
        p_out[13] = phi_perturb
        
        # Parameter names
        pnames_out[11] = 'A_perturb'
        pnames_out[12] = 'theta_perturb'
        pnames_out[13] = 'phi_perturb'
    elif PTC_or_ISO == 'ISO':
        # Parameter vector
        p_out[11] = d_x
        p_out[12] = d_y
        p_out[13] = d_z
        p_out[14] = iso1
        p_out[15] = iso2
        p_out[16] = iso3
        
        # Parameter names
        pnames_out[11] = 'd_x'
        pnames_out[12] = 'd_y'
        pnames_out[13] = 'd_z'
        pnames_out[14] = 'iso1'
        pnames_out[15] = 'iso2'
        pnames_out[16] = 'iso3'
        
    return x_init_out, p_out, pnames_out