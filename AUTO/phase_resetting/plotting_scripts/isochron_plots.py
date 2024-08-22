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
    from plotting_scripts.initial_PO_plots import plot_base_periodic_orbit
    from scipy.io import loadmat
    
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
    iso3 = run_in['iso3']

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Single Isochron Test', figsize=[8, 8])
    ax = plt.axes(projection='3d')

    #--------------#
    #     Plot     #
    #--------------#
    # Plot phase space portrait of periodic orbit
    ax = plot_base_periodic_orbit(ax)

    # Plot isochron
    ax.plot(iso1, iso2, iso3, color='C0')
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # ax.set_xticks(arange(0, 12, 2))
    # ax.set_xticks(arange(1, 12, 2), minor=True)
    
    # ax.set_yticks(arange(0, 10, 2))
    # ax.set_yticks(arange(1, 10, 2), minor=True)

    # ax.set_zticks(arange(0.0, 22, 4.0))
    # ax.set_zticks(arange(2.0, 22, 4.0), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(0.0, 8.0)
    ax.set_zlim(0.0, 20)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$G(t)$')
    ax.set_ylabel(r'$Q(t)$')
    ax.set_zlabel(r'$I(t)$')
    
    ax.set_title(r'Single Isochron (TEST)')
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid(False)
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    fig.savefig(filename_out)
    fig.show()