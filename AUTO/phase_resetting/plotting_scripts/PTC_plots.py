#------------------------------------------------------------------------------#
# Plot single isochron against the periodic orbit in phase space
def plot_single_PTC(run_in):
    """
    Plots a single PTC from the continuation run "run08_PTC_single"

    Input
    -------
    run_in : AUTO generated run file
        The continuation run (probably run08_isochron_single).
    """
    import matplotlib.pyplot as plt
    from numpy import arange, cos, sin
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # Filename for figure
    filename_out = "./images/PTC_single.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Cooridnates of the isochron
    theta_old = run_in['theta_old']
    theta_new = run_in['theta_new']

    # Solution to read
    sol_read = run_in(1)

    # Perturbation amplitude
    A_perturb     =  sol_read['A_perturb']
    # Perturbation direction angle
    theta_perturb = sol_read['theta_perturb']
    # Perturbation vector
    d_perturb     = [cos(theta_perturb), sin(theta_perturb)]

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Single PTC Test', figsize=[8, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot isochron
    ax.plot(theta_old, theta_new)
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(0.0, 2.5, 0.5))
    ax.set_xticks(arange(0.25, 2.5, 0.5), minor=True)
    
    ax.set_yticks(arange(0.0, 2.5, 0.5))
    ax.set_yticks(arange(0.25, 2.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 2.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$\theta_{\mathrm{old}}$')
    ax.set_ylabel(r'$\theta_{\mathrm{new}}$')
    ax.set_title((r'Single PTC (TEST) with $A_{{\mathrm{{p}}}} = {:.2f}$ and '
                  r'$\vec{{d}}_{{\mathrm{{p}}}} = \left( {:.1f}, {:.1f} \right)$'
                  ).format(A_perturb, d_perturb[0], d_perturb[1]))
    
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

#------------------------------------------------------------------------------#
def read_PTC_scan_data(run_str_in):
    from python_files.write_data_functions import bd_read

    # Get sorted data folders
    data_dir = sort_data_folders(run_str_in)

    # Empty arrays for data
    theta_old = []
    theta_new = []
    A_perturb = []

    # Read data
    for idx, run in enumerate(data_dir):
        # Read data
        bd = bd_read('{}/{}'.format(run_str_in, run))

        # Read iso1 values
        theta_old_read = bd['theta_old']

        # Read iso1 values
        theta_new_read = bd['theta_new']

        # Read A_perturb value
        sol = bd(1)
        A_perturb_read = sol['A_perturb']

        # Append data
        theta_old.append(theta_old_read)
        theta_new.append(theta_new_read)
        A_perturb.append(A_perturb_read)

    #----------------#
    #     Output     #
    #----------------#
    return theta_old, theta_new, A_perturb

#------------------------------------------------------------------------------#
# Plot single isochron against the periodic orbit in phase space
def plot_multi_PTC(run_str_in):
    """
    Plots all isochrons from the continuation scan run "run09_isochron_multi"

    Input
    -------
    run_str_in : str
        String label identifier for the current run.
    """
    import matplotlib.pyplot as plt
    from python_files.write_data_functions import bd_read
    from numpy import arange, cos, sin
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # Filename for figure
    filename_out = "./images/PTC_multi.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Cooridnates of the isochron
    theta_old, theta_new, A_perturb = read_PTC_scan_data(run_str_in)

    run_in = bd_read('{}/sol_001'.format(run_str_in))
    sol_read = run_in(1)

    # Solution to read
    sol_read = run_in(1)

    # Perturbation direction angle
    theta_perturb = sol_read['theta_perturb']
    # Perturbation vector
    d_perturb     = [cos(theta_perturb), sin(theta_perturb)]

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Multiple PTCs', figsize=[8, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot PTC
    for i in range(len(theta_old)):
        old_plot = theta_old[i]
        new_plot = theta_new[i]

        # Set colour and linestlye
        if i <= 9:
            colour = 'C{}'.format(i)
            ls     = 'solid'
        elif i > 9 and i <= 19:
            colour = 'C{}'.format(i-10)
            ls     = 'dashed'
        elif i > 19 and i <= 29:
            colour = 'C{}'.format(i-20)
            ls     = 'dotted'

        # Plot without label
        ax.plot(old_plot, new_plot, color=colour, ls=ls,
                label=r'$A_{{\mathrm{{p}}}} = {:.3f}$'.format(A_perturb[i]))
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(0.0, 2.5, 0.5))
    ax.set_xticks(arange(0.25, 2.5, 0.5), minor=True)
    
    ax.set_yticks(arange(0.0, 2.5, 0.5))
    ax.set_yticks(arange(0.25, 2.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 2.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$\theta_{\mathrm{old}}$')
    ax.set_ylabel(r'$\theta_{\mathrm{new}}$')
    
    ax.set_title((r'PTC with'
                  r'$\vec{{d}}_{{\mathrm{{p}}}} = \left( {:.1f}, 0.0, {:.1f} \right)$'
                  ).format(d_perturb[0], d_perturb[1]))
    
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
# Plot single level set (A_perturb against theta_perturb)
def plot_single_level_set(run_in):
    """
    Plots a single level set of A_perturb vs theta_perturb for a fixed
    theta_new and theta_old.

    Input
    -------
    run_in : AUTO generated run file
        The continuation run (probably run08_isochron_single).
    """
    import matplotlib.pyplot as plt
    from numpy import arange, cos, sin, pi
    plt.close('all')
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # Filename for figure
    filename_out = "./images/PTC_single.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Perturbation amplitude
    A_perturb     =  run_in['A_perturb']
    # Perturbation direction angle
    theta_perturb = run_in['theta_perturb']

    # Solution to read
    sol_read = run_in(1)

    # theta_old and theta_new
    theta_old = sol_read['theta_old']
    theta_new = sol_read['theta_new']

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Single Level Set (Test)', figsize=[8, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot isochron
    ax.plot(theta_perturb / (pi), A_perturb)
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # ax.set_xticks(arange(0.0, 2.5, 0.5))
    # ax.set_xticks(arange(0.25, 2.5, 0.5), minor=True)
    
    # ax.set_yticks(arange(0.0, 2.5, 0.5))
    # ax.set_yticks(arange(0.25, 2.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    # ax.set_xlim(0.0, 2.0)
    # ax.set_ylim(0.0, 2.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$\theta_{\mathrm{p}} / \pi$')
    ax.set_ylabel(r'$A_{\mathrm{p}}$')

    # Axis title
    ax.set_title((r'Single Level Set (TEST) with $\theta_{{\mathrm{{old}}}} = {:.2f}$, '
                  r'$\theta_{{\mathrm{{new}}}} = {:.2f}$'
                  ).format(theta_old, theta_new))
    
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
def read_level_set_scan_data(run_str_in):
    from python_files.write_data_functions import bd_read

    # Get sorted data folders
    data_dir = sort_data_folders(run_str_in)

    # Empty arrays for data
    theta_old = []
    theta_new = []
    A_perturb = []

    # Read data
    for idx, run in enumerate(data_dir):
        # Read data
        bd = bd_read('{}/{}'.format(run_str_in, run))

        # Read iso1 values
        theta_old_read = bd['theta_old']

        # Read iso1 values
        theta_new_read = bd['theta_new']

        # Read A_perturb value
        sol = bd(1)
        A_perturb_read = sol['A_perturb']

        # Append data
        theta_old.append(theta_old_read)
        theta_new.append(theta_new_read)
        A_perturb.append(A_perturb_read)

    #----------------#
    #     Output     #
    #----------------#
    return theta_old, theta_new, A_perturb

#------------------------------------------------------------------------------#
# Plot single isochron against the periodic orbit in phase space
def plot_multi_level_set(run_str_in):
    """
    Plots all level sets from the continuation scan run "run09_isochron_multi"

    Input
    -------
    run_str_in : str
        String label identifier for the current run.
    """
    import matplotlib.pyplot as plt
    from python_files.write_data_functions import bd_read
    from numpy import arange, cos, sin
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # Filename for figure
    filename_out = "./images/PTC_multi.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Cooridnates of the isochron
    theta_old, theta_new, A_perturb = read_PTC_scan_data(run_str_in)

    run_in = bd_read('{}/sol_001'.format(run_str_in))
    sol_read = run_in(1)

    # Solution to read
    sol_read = run_in(1)

    # Perturbation direction angle
    theta_perturb = sol_read['theta_perturb']
    # Perturbation vector
    d_perturb     = [cos(theta_perturb), sin(theta_perturb)]

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Multiple PTCs', figsize=[8, 8])
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot PTC
    for i in range(len(theta_old)):
        old_plot = theta_old[i]
        new_plot = theta_new[i]

        # Set colour and linestlye
        if i <= 9:
            colour = 'C{}'.format(i)
            ls     = 'solid'
        elif i > 9 and i <= 19:
            colour = 'C{}'.format(i-10)
            ls     = 'dashed'
        elif i > 19 and i <= 29:
            colour = 'C{}'.format(i-20)
            ls     = 'dotted'

        # Plot without label
        ax.plot(old_plot, new_plot, color=colour, ls=ls,
                label=r'$A_{{\mathrm{{p}}}} = {:.3f}$'.format(A_perturb[i]))
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(0.0, 2.5, 0.5))
    ax.set_xticks(arange(0.25, 2.5, 0.5), minor=True)
    
    ax.set_yticks(arange(0.0, 2.5, 0.5))
    ax.set_yticks(arange(0.25, 2.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 2.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$\theta_{\mathrm{old}}$')
    ax.set_ylabel(r'$\theta_{\mathrm{new}}$')
    
    ax.set_title((r'PTC with'
                  r'$\vec{{d}}_{{\mathrm{{p}}}} = \left( {:.1f}, 0.0, {:.1f} \right)$'
                  ).format(d_perturb[0], d_perturb[1]))
    
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
