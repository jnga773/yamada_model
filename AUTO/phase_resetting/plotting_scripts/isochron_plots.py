#------------------------------------------------------------------------------#
# Plot single isochron against the periodic orbit in phase space
def plot_single_isochron(run_in):
    """
    Plots a single isochron from the continuation run "run08_isochron_single"

    Input
    -------
    run_in : AUTO generated run file
        The continuation run (probably run08_isochron_single).
    """
    import matplotlib.pyplot as plt
    from numpy import arange
    from plotting_scripts.initial_plots import plot_phase_space_template
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # Filename for figure
    filename_out = "./images/isochron_single.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Cooridnates of the isochron
    iso1 = run_in['iso1']
    iso2 = run_in['iso2']

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Single Isochron Test', figsize=[8, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot phase space portrait of periodic orbit
    ax = plot_phase_space_template(ax)

    # Plot isochron
    ax.plot(iso1, iso2, color='C0')
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(-2.0, 2.5, 0.5))
    ax.set_xticks(arange(-1.75, 2.5, 0.5), minor=True)
    
    ax.set_yticks(arange(-2.0, 2.5, 0.5))
    ax.set_yticks(arange(-1.75, 2.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-2.0, 2.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$x_{1}(t)$')
    ax.set_ylabel(r'$x_{2}(t)$')
    ax.set_title(r'Single Isochron (TEST)')
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid()
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    fig.savefig(filename_out)
    fig.show()

#------------------------------------------------------------------------------#
def sort_data_folders(run_str_in):
    from os import listdir
    from numpy import argsort
    
    # List all directories in ./data/run08_isochron_scan/
    dir_main = './data/{}/'.format(run_str_in)
    dir_list = listdir(dir_main)

    dir_list = [f for f in dir_list if not f.startswith('.')]
    
    # Cyycle through and get the integer number
    int_list = []
    for f in dir_list:
        int_read = f[4:]
        int_list.append(int(int_read))
        
    # Sort int_list and get indices
    idx_sorted = argsort(int_list)
    
    dir_list_sorted = [dir_list[i] for i in idx_sorted]
    
    return dir_list_sorted

def read_isochron_scan_data(run_str_in):
    from python_files.write_data_functions import bd_read

    # Get sorted data folders
    data_dir = sort_data_folders(run_str_in)

    # Empty arrays for data
    iso1 = []
    iso2 = []

    # Read data
    for idx, run in enumerate(data_dir):
        # Read data
        bd = bd_read('{}/{}'.format(run_str_in, run))

        # Read iso1 values
        iso1_read = bd['iso1']

        # Read iso1 values
        iso2_read = bd['iso2']

        # Append data
        iso1.append(iso1_read)
        iso2.append(iso2_read)

    #----------------#
    #     Output     #
    #----------------#
    return iso1, iso2

# Plot single isochron against the periodic orbit in phase space
def plot_multi_isochron(run_str_in):
    """
    Plots all isochrons from the continuation scan run "run09_isochron_multi"

    Input
    -------
    run_str_in : str
        String label identifier for the current run.
    """
    import matplotlib.pyplot as plt
    from numpy import arange
    from plotting_scripts.initial_plots import plot_phase_space_template
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # Filename for figure
    filename_out = "./images/isochron_multi.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Cooridnates of the isochron
    iso1, iso2 = read_isochron_scan_data(run_str_in)

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Multplie Isochrons', figsize=[8, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot phase space portrait of periodic orbit
    ax = plot_phase_space_template(ax)

    # Plot isochrons
    for i in range(len(iso1)):
        iso1_plot = iso1[i]
        iso2_plot = iso2[i]

        if i == 0:
            # Plot with a label
            ax.plot(iso1_plot, iso2_plot, color='C0', alpha=0.5,
                    label='isochron')
        else:
            # Plot without label
            ax.plot(iso1_plot, iso2_plot, color='C0', alpha=0.5)
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(-1.5, 2.5, 0.5))
    ax.set_xticks(arange(-1.25, 2.5, 0.5), minor=True)
    
    ax.set_yticks(arange(-1.5, 2.5, 0.5))
    ax.set_yticks(arange(-1.25, 2.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(-1.5, 2.0)
    ax.set_ylim(-1.5, 2.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$x_{1}(t)$')
    ax.set_ylabel(r'$x_{2}(t)$')
    ax.set_title(r'Multiple Isochrons')
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid()
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    fig.savefig(filename_out)
    fig.show()
