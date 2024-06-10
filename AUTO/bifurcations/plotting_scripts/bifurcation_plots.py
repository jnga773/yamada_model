#------------------------------------------------------------------------------#
# Data to read bifurcation data and write to MATLAB .mat file
def save_data_mat(run_names_in):
    """
    Reads the A and \gamma values from the various bifurcation data runs,
    then saves it to a MATLAB .mat file.

    Input
    -------
    run_names_in : dict
      Dictionary of all of the different string identifiers for the bifurcation
      data to be read. Must contain the following keys:
          - 'H' : Hopf bifurcation data,
          - 'S' : Saddle bifurcation,
          - 'T' : Transcritical bifurcations,
          - 'L' : Homoclinic,
          - 'D' : Double-limit cycles.
    """
    from python_files.write_data_functions import bd_read
    from numpy import argmin
    from scipy.io import savemat    

    #------------------------------#
    #     Read Data: H and NSA     #
    #------------------------------#
    # Read bifurcation data
    bd = bd_read(run_names_in['H'])
    
    # Read data
    A_read     = bd['A']
    gamma_read = bd['gamma']

    # Cut off gamma < 0
    mask       = gamma_read >= 0.0
    A_read     = A_read[mask]
    gamma_read = gamma_read[mask]

    # Find point where A = min(A)
    min_idx = argmin(A_read)

    # Split into Hopf and neutral saddles
    A_SQ     = A_read[:min_idx]
    gamma_SQ = gamma_read[:min_idx]
    
    A_H     = A_read[min_idx:]
    gamma_H = gamma_read[min_idx:]

    #------------------------#
    #     Read Data: A_S     #
    #------------------------#
    # Read bifurcation data
    bd = bd_read(run_names_in['S'])

    # Read data
    A_S     = bd['A']
    gamma_S = bd['gamma']

    #------------------------#
    #     Read Data: A_T     #
    #------------------------#
    # Read bifurcation data
    bd = bd_read(run_names_in['T'])

    # Read data
    A_T     = bd['A']
    gamma_T = bd['gamma']

    #----------------------#
    #     Read Data: D     #
    #----------------------#
    # Read bifurcation data
    bd = bd_read(run_names_in['D'])

    # Read data
    A_D     = bd['A']
    gamma_D = bd['gamma']

    #----------------------#
    #     Read Data: L     #
    #----------------------#
    # Read bifurcation data
    bd = bd_read(run_names_in['L'])

    # Read data
    A_L     = bd['A']
    gamma_L = bd['gamma']

    #-------------------#
    #     Save Data     #
    #-------------------#
    # Output dictionary
    dict_out = {'A_H' : A_H , 'gamma_H' : gamma_H,
                'A_SQ': A_SQ, 'gamma_SQ': gamma_SQ,
                'A_S' : A_S , 'gamma_S' : gamma_S,
                'A_T' : A_T , 'gamma_T' : gamma_T,
                'A_D' : A_D , 'gamma_D' : gamma_D,
                'A_L' : A_L , 'gamma_L' : gamma_L}
    
    # Save matrix file
    savemat('./plotting_scripts/bifurcation_data.mat', dict_out)

#------------------------------------------------------------------------------#
def read_mat_data(data_in, file_in):
    # Read data
    A_read     = data_in['A_{}'.format(file_in)]
    gamma_read = data_in['gamma_{}'.format(file_in)]

    # Reshape
    A_out     = A_read.reshape(A_read.shape[1])
    gamma_out = gamma_read.reshape(gamma_read.shape[1])

    return A_out, gamma_out

#------------------------------------------------------------------------------#
# Plot regular bifurcation diagram
def plot_bifurcation_diagram():
    """
    Plots the two-paramater (A, \gamma) bifurcation diagram of the Yamada model.
    """
    from scipy.io import loadmat
    import matplotlib.pyplot as plt
    from numpy import arange

    # Add style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    #-------------------#
    #     Read Data     #
    #-------------------#
    data = loadmat('./plotting_scripts/bifurcation_data.mat')

    # Hopf bifurcation data
    A_H, gamma_H = read_mat_data(data, 'H')

    # Neutral saddle=node
    A_SQ, gamma_SQ = read_mat_data(data, 'SQ')

    # Saddle-Node data
    A_S, gamma_S = read_mat_data(data, 'S')

    # Transcritical data
    A_T, gamma_T = read_mat_data(data, 'T')

    # Homoclinic data
    A_L, gamma_L = read_mat_data(data, 'L')

    # Double limit cycle
    A_D, gamma_D = read_mat_data(data, 'D')

    #--------------------------------------------------------------------------#
    #                            PLOT (Full Picture)                           #
    #--------------------------------------------------------------------------#
    fig = plt.figure(num='Bifurcation Diagram', figsize=[12, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    ax.plot(A_S, gamma_S, color='C0', ls='solid', label='S')
    ax.plot(A_T, gamma_T, color='C1', ls='dashed', label='T')
    ax.plot(A_H, gamma_H, color='C2', ls='solid', label='H')
    ax.plot(A_SQ, gamma_SQ, color='C3', ls='solid', label='SQ')
    ax.plot(A_L, gamma_L, color='C4', ls='solid', label='L')
    ax.plot(A_D, gamma_D, color='C5', ls='solid', label='D')

    # Legend
    ax.legend()

    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # Ticks - X axis
    ax.set_xticks(arange(5.0, 12.0, 1.0))
    ax.set_xticks(arange(5.2, 11.2, 0.2), minor=True)

    # Ticks - Y Axis
    ax.set_yticks(arange(0.0, 0.30, 0.05))
    ax.set_yticks(arange(0.01, 0.25, 0.01), minor=True)

    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(5.0, 11.0)
    ax.set_ylim(0.0, 0.25)

    #--------------#
    #     Grid     #
    #--------------#
    ax.grid()

    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$A$')
    ax.set_ylabel(r'$\gamma$')

    # Title
    ax.set_title(r'Yamada Model with $\left(B = 5.80, a = 1.80 \right)$')

    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout(pad=0)
    fig.savefig("./images/bifurcations.pdf")
    fig.show()

#------------------------------------------------------------------------------#
# Plot zoomed bifurcation diagram
def plot_bifurcation_diagram_zoomed():
    """
    Plots the two-paramater (A, \gamma) bifurcation diagram of the Yamada model.
    """
    from scipy.io import loadmat
    import matplotlib.pyplot as plt
    from numpy import arange

    # Add style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    #-------------------#
    #     Read Data     #
    #-------------------#
    data = loadmat('./plotting_scripts/bifurcation_data.mat')

    # Hopf bifurcation data
    A_H, gamma_H = read_mat_data(data, 'H')

    # Neutral saddle=node
    A_SQ, gamma_SQ = read_mat_data(data, 'SQ')

    # Saddle-Node data
    A_S, gamma_S = read_mat_data(data, 'S')

    # Transcritical data
    A_T, gamma_T = read_mat_data(data, 'T')

    # Homoclinic data
    A_L, gamma_L = read_mat_data(data, 'L')

    # Double limit cycle
    A_D, gamma_D = read_mat_data(data, 'D')

    #--------------------------------------------------------------------------#
    #                           PLOT (Zoomed Picture)                          #
    #--------------------------------------------------------------------------#
    fig = plt.figure(num='Bifurcation Diagram (Zoomed)', figsize=[12, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    ax.plot(A_S, gamma_S, color='C0', ls='solid', label='S')
    ax.plot(A_T, gamma_T, color='C1', ls='dashed', label='T')
    ax.plot(A_H, gamma_H, color='C2', ls='solid', label='H')
    ax.plot(A_SQ, gamma_SQ, color='C3', ls='solid', label='SQ')
    ax.plot(A_L, gamma_L, color='C4', ls='solid', label='L')
    ax.plot(A_D, gamma_D, color='C5', ls='solid', label='D')

    # Legend
    ax.legend()

    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # Ticks - X axis
    ax.set_xticks(arange(6.6, 7.6, 0.2))
    ax.set_xticks(arange(6.5, 7.45, 0.05), minor=True)

    # Ticks - Y Axis
    ax.set_yticks(arange(0.04, 0.11, 0.01))
    ax.set_yticks(arange(0.042, 0.12, 0.002), minor=True)

    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(6.5, 7.4)
    ax.set_ylim(0.04, 0.10)

    #--------------#
    #     Grid     #
    #--------------#
    ax.grid()

    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$A$')
    ax.set_ylabel(r'$\gamma$')

    # Title
    ax.set_title(r'Yamada Model with $\left(B = 5.80, a = 1.80 \right)$')

    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout(pad=0)
    fig.savefig("./images/bifurcations_zoomed.pdf")
    fig.show()