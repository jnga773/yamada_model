#==============================================================================#
#                FUNCTIONS USED IN THE PHASE RESET CALCULATIONS                #
#==============================================================================#
# This Python module contains the following functions:
#
# - calc_initial_solution_phase_reset : Calculates and sets the initial
#                                       solution to solve for the phase
#                                       resetting problems, and write it
#                                       to ./data/[XXX].dat file.

#------------------------------------------------------------------------------#
# Calculate the initial solutions to the phase resetting problem.
def write_initial_solution_phase_reset(k_in):
    """
    Calculates and sets the initial solution to solve for the
    phase resetting problem and write it to a .dat file
    """
    from scipy.io import loadmat

    #---------------------------------------------------------------------------#
    #                                 Read Data                                 #
    #---------------------------------------------------------------------------#
    # Load data from .mat file
    data_floquet = loadmat('./data_mat/floquet_solution.mat')

    # Read time data
    t_read = data_floquet['t_read'].flatten()

    #--------------------#
    #     State data     #
    #--------------------#
    # Periodic orbit solution
    x1_read = data_floquet['x1_read'].flatten()
    x2_read = data_floquet['x2_read'].flatten()
    x3_read = data_floquet['x3_read'].flatten()

    # Initial point
    x1_0 = x1_read[0]
    x2_0 = x2_read[0]
    x3_0 = x3_read[0]

    #------------------------------#
    #     Perpindicular vector     #
    #------------------------------#
    wn1_read = data_floquet['wn1_read'].flatten()
    wn2_read = data_floquet['wn2_read'].flatten()
    wn3_read = data_floquet['wn3_read'].flatten()

    # Initial point
    wn1_0 = wn1_read[0]
    wn2_0 = wn2_read[0]
    wn3_0 = wn3_read[0]

    #---------------------------------------------------------------------------#
    #                              "Periodise" Data                             #
    #---------------------------------------------------------------------------#
    from numpy import concatenate

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

    #---------------------------------------------------------------------------#
    #                              Interpolate Data                             #
    #---------------------------------------------------------------------------#
    from numpy import interp, ones

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
            
    #---------------------------------------------------------------------------#
    #                                Write Data                                 #
    #---------------------------------------------------------------------------#
    from numpy import array
    
    # Time
    t_seg = t_period

    # Segment 1
    x_seg1 = array([x1_interp, x2_interp, x3_interp,
                    wn1_interp, wn2_interp, wn3_interp])
    x_seg1 = x_seg1.T

    # Segment 2
    x_seg2 = array([x1_extended, x2_extended, x3_extended,
                    wn1_extended, wn2_extended, wn3_extended])
    x_seg2 = x_seg2.T

    # Segment 3
    x_seg3 = array([x1_extended, x2_extended, x3_extended])
    x_seg3 = x_seg3.T

    # Segment 4
    x_seg4 = array([x1_period, x2_period, x3_period])
    x_seg4 = x_seg4.T

    filename_data = './data_mat/initial_solution_PR.dat'
    with open(filename_data, 'w') as file_data:
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

#------------------------------------------------------------------------------#