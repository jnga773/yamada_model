#------------------------------------------------------------------------------#
# Plot initial phase space template
def plot_base_periodic_orbit(ax_in, plot_dimension='3D'):
    """
    Plots the initial phase space portrait of the periodic orbit from the saved
    data in ./data/.

    Input
    -------
    ax_in : matplotlib.pyplot ax structure
        Input axis to plot the solutions.
    2D_or_3D : str
        Plots the 3D phase space portrait, or 2D in the (G-I) plane.

    Output
    -------
    ax_in : The same axis things but with plots on it now
    """
    from scipy.io import loadmat
    import matplotlib.pyplot as plt

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Load data from matrix
    data = loadmat('./data_mat/initial_PO.mat')

    # Plot things
    xbp  = data['xbp_PO']
    # x0   = data['x0']
    xpos = data['xpos']
    # xneg = data['xneg']

    #--------------#
    #     Plot     #
    #--------------#
    if plot_dimension == '3D':
        # Plot periodic orbit
        ax_in.plot(xbp[:, 0], xbp[:, 1], xbp[:, 2], color='k', ls='solid', label=r'$\Gamma$')
        
        # # Plot x0
        # ax_in.plot(x0[0], x0[1], x0[2], color='r', ls='none', label=r'$o$',
        #           marker='o', markeredgecolor='r', markerfacecolor='r',
        #           markersize=12)
        
        # Plot xpos
        ax_in.plot(xpos[0], xpos[1], xpos[2], color='b', ls='none', label=r'$q$',
                marker='*', markeredgecolor='b', markerfacecolor='b',
                markersize=12)
        
        # # Plot xneg
        # ax_in.plot(xneg[0], xneg[1], xneg[2], color='r', ls='none', label=r'$p$',
        #           marker='*', markeredgecolor='r', markerfacecolor='r',
        #           markersize=12)

    elif plot_dimension == '2D':
        # Plot periodic orbit
        ax_in.plot(xbp[:, 0], xbp[:, 2], color='k', ls='solid', label=r'$\Gamma$')
        
        # Plot xpos
        ax_in.plot(xpos[0], xpos[2], color='b', ls='none', label=r'$q$',
                marker='*', markeredgecolor='b', markerfacecolor='b',
                markersize=12)

#------------------------------------------------------------------------------#
# Plotting function template for phase resetting problem for each segment
def plot_PR_solution(ax_in, sol_in, segment, plot_dimension='3D'):
    """
    Plots the base solution for each segment of the phase-resetting problem.
    """
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read state-space solution
    if segment == 'seg1':
        # Read solution
        x1 = sol_in['{}_x1'.format(segment)]
        x2 = sol_in['{}_x2'.format(segment)]
        x3 = sol_in['{}_x3'.format(segment)]

        # Plot colour
        colour = 'C0'

    elif segment == 'seg2':
        # Read solution
        x1 = sol_in['{}_x1'.format(segment)]
        x2 = sol_in['{}_x2'.format(segment)]
        x3 = sol_in['{}_x3'.format(segment)]

        # Plot colour
        colour = 'C1'

    elif segment == 'seg3':
        # Read solution
        x1 = sol_in['{}_x1'.format(segment)]
        x2 = sol_in['{}_x2'.format(segment)]
        x3 = sol_in['{}_x3'.format(segment)]

        # Plot colour
        colour = 'C2'

    elif segment == 'seg4':
        # Read solution
        x1 = sol_in['{}_x1'.format(segment)]
        x2 = sol_in['{}_x2'.format(segment)]
        x3 = sol_in['{}_x3'.format(segment)]

        # Plot colour
        colour = 'C3'

    #--------------#
    #     Plot     #
    #--------------#
    # Plot the state-space solution
    if plot_dimension == '3D':
        # Plot state solutions
        ax_in.plot(x1, x2, x3, color=colour, alpha=0.5, label=segment)
        # Plot initial point
        ax_in.plot(x1[0], x2[0], x3[0], ls='none', marker='s',
                   markerfacecolor=colour, markeredgecolor=colour,
                   label=rf'$x_{{{segment[-1]}}}(0)$')
        # Plot final point
        ax_in.plot(x1[-1], x2[-1], x3[-1], ls='none', marker='d',
                   markerfacecolor=colour, markeredgecolor=colour,
                   label=rf'$x_{{{segment[-1]}}}(1)$')

    elif plot_dimension == '2D':
        # Plot state solutions
        ax_in.plot(x1, x3, color=colour, alpha=0.5, label=segment)
        # Plot initial point
        ax_in.plot(x1[0], x3[0], ls='none', marker='s',
                   markerfacecolor=colour, markeredgecolor=colour,
                   label=rf'$x_{{{segment[-1]}}}(0)$')
        # Plot final point
        ax_in.plot(x1[-1], x3[-1], ls='none', marker='d',
                   markerfacecolor=colour, markeredgecolor=colour,
                   label=rf'$x_{{{segment[-1]}}}(1)$')

#------------------------------------------------------------------------------#
# Plotting function for phase resetting problem
def plot_phase_reset_phase_space(run_in, label_in):
    """
    Plots the phase-space solution of the shifted PO BVP continuation run from
    the saved MATLAB .mat data file.
    """
    from scipy.io import loadmat
    import matplotlib.pyplot as plt
    from numpy import arange
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')
    
    # Figure name for file
    filename_out = './images/phase_reset_curve.pdf'
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    # Load data from matrix
    data = loadmat('./data_mat/initial_PO.mat')

    # Plot things
    xbp = data['xbp_PO']
    x0  = data['x0']

    # Phase reset solution
    sol = run_in(label_in)

    # PTC parameters
    theta_old = sol['theta_old']
    theta_new = sol['theta_new']
    A_perturb = sol['A_perturb']
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    fignum = 'Phase Reset Periodic Orbit (3D)'
    plt.close(fignum)
    fig = plt.figure(num=fignum, figsize=[10, 10])
    ax = plt.axes(projection='3d')
    
    #--------------#
    #     Plot     #
    #--------------#
    # Plot unperturbed periodic orbit
    plot_base_periodic_orbit(ax, '3D')

    # Plot phase-resetting segments
    # plot_PR_solution(ax, sol, 'seg1', '3D')
    # plot_PR_solution(ax, sol, 'seg2', '3D')
    # plot_PR_solution(ax, sol, 'seg3', '3D')
    plot_PR_solution(ax, sol, 'seg4', '3D')
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(0, 12, 2))
    ax.set_xticks(arange(1, 12, 2), minor=True)
    
    ax.set_yticks(arange(0, 10, 2))
    ax.set_yticks(arange(1, 10, 2), minor=True)

    ax.set_zticks(arange(0.0, 22, 4.0))
    ax.set_zticks(arange(2.0, 22, 4.0), minor=True)
    
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

    ax.set_title((r'$\theta_{{\mathrm{{old}}}} = {:.3f}, \theta_{{\mathrm{{new}}}} = {:.3f}, A_{{\mathrm{{perturb}}}} = {:.3f}$'
                  ).format(theta_old, theta_new, A_perturb))
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid(which='major')
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    fig.savefig(filename_out)
    fig.show()

#------------------------------------------------------------------------------#
# Plotting function for phase resetting problem
def plot_phase_reset_phase_space_2D(run_in, label_in):
    """
    Plots the phase-space solution of the shifted PO BVP continuation run from
    the saved MATLAB .mat data file.
    """
    from scipy.io import loadmat
    import matplotlib.pyplot as plt
    from numpy import arange
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')
    
    # Figure name for file
    filename_out = './images/phase_reset_curve.pdf'
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    # Load data from matrix
    data = loadmat('./data_mat/initial_PO.mat')

    # Plot things
    xbp = data['xbp_PO']
    x0  = data['x0']

    # Phase reset solution
    sol = run_in(label_in)
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    fignum = 'Phase Reset Periodic Orbit (2D)'
    plt.close(fignum)
    fig = plt.figure(num=fignum, figsize=[10, 10])
    ax = plt.gca()
    
    #--------------#
    #     Plot     #
    #--------------#
    # Plot unperturbed periodic orbit
    # plot_base_periodic_orbit(ax, '2D')

    # Plot phase-resetting segments
    # plot_PR_solution(ax, sol, 'seg1', '2D')
    # plot_PR_solution(ax, sol, 'seg2', '2D')
    # plot_PR_solution(ax, sol, 'seg3', '2D')
    plot_PR_solution(ax, sol, 'seg4', '2D')
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(0, 12, 2))
    ax.set_xticks(arange(1, 12, 2), minor=True)

    ax.set_yticks(arange(0.0, 22, 4.0))
    ax.set_yticks(arange(2.0, 22, 4.0), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(0.0, 20)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$G(t)$')
    ax.set_ylabel(r'$I(t)$')

    ax.set_title(r'Shifted Periodic Orbit')
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid(which='major')
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    fig.savefig(filename_out)
    fig.show()

#------------------------------------------------------------------------------#
# Plot the perturbed orbit segment (seg4)
def plot_perturbed_states(run_in, label_in):
    """
    Plots the perturbed states against the unperturbed orbit.
    """
    from scipy.io import loadmat
    import matplotlib.pyplot as plt
    from numpy import arange, mod, concatenate
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    #--------------------------------------------------------------------------#
    #                                Read Data                                 #
    #--------------------------------------------------------------------------#
    #-------------------------------#
    #     Data: Perturbed Orbit     #
    #-------------------------------#
    # Read phase-reset data
    sol = run_in(label_in)

    # Periodicity
    k = int(sol['k'])

    # State space solution
    xbp = [sol['seg4_x1'], sol['seg4_x2'], sol['seg4_x3']]

    # Time solution
    tbp = sol['t'] * k

    # Add \theta_old to time array for point where perturbation is applied
    theta_old = sol['theta_old']
    if theta_old < 1.0:
        tbp += theta_old
    else:
        tbp += mod(theta_old, 1.0)

    # PTC parameters
    theta_old = sol['theta_old']
    theta_new = sol['theta_new']
    A_perturb = sol['A_perturb']

    #---------------------------------#
    #     Data: Unperturbed Orbit     #
    #---------------------------------#
    # Load in the data of the unperturbed orbit
    data_PO = loadmat('./data_mat/initial_PO.mat')

    # State-space solution
    xbp_PO  = data_PO['xbp_PO']
    # Time data
    tbp_PO  = data_PO['tbp'].flatten()
    tbp_PO *= 1 / max(tbp_PO)

    # Append k periods of data
    x1_plot  = xbp_PO[:, 0]
    x2_plot  = xbp_PO[:, 1]
    x3_plot  = xbp_PO[:, 2]
    tbp_plot = tbp_PO
    t_max    = tbp_plot[-1]

    for i in range(k-1):
        # Append time
        tbp_plot = concatenate((tbp_plot, t_max + tbp_PO[1:]))
        t_max = tbp_plot[-1]

        # Append states
        x1_plot = concatenate((x1_plot, xbp_PO[1:, 0]))
        x2_plot = concatenate((x2_plot, xbp_PO[1:, 1]))
        x3_plot = concatenate((x3_plot, xbp_PO[1:, 2]))

    # Merge xbp data
    xbp_plot = [x1_plot, x2_plot, x3_plot]

    #--------------------------------------------------------------------------#
    #                                Plot Data                                 #
    #--------------------------------------------------------------------------#
    fignum = 'Perturbed States'
    plt.close(fignum)
    fig, ax = plt.subplots(ncols=1, nrows=3, sharex=True,
                           num=fignum, figsize=[12, 8])
    
    #--------------#
    #     Plot     #
    #--------------#
    # Plot
    for i in range(3):
        ax[i].plot(tbp_plot, xbp_plot[i], color='C1', ls='solid', lw=2.0, alpha=0.5)
        ax[i].plot(tbp, xbp[i], color='C0', ls='solid', lw=1.0)

    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax[2].set_xticks(arange(0.0, k+5, 5))
    ax[2].set_xticks(arange(2.5, k+5, 5), minor=True)

    ax[0].set_yticks(arange(0.0, 6.0, 1.0))
    ax[0].set_yticks(arange(0.5, 6.0, 1.0), minor=True)

    ax[1].set_yticks(arange(0.0, 5.0, 1.0))
    ax[1].set_yticks(arange(0.5, 5.0, 1.0), minor=True)

    ax[2].set_yticks(arange(0.0, 30.0, 5.0))
    ax[2].set_yticks(arange(2.5, 30.0, 5.0), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax[2].set_xlim(-0.01, k+0.1)
    
    ax[0].set_ylim(0.0, 5.0)
    ax[1].set_ylim(0.0, 4.0)
    ax[2].set_ylim(0.0, 25.0)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax[2].set_xlabel(r'$t$')

    ax[0].set_ylabel(r'$G(t)$')
    ax[1].set_ylabel(r'$Q(t)$')
    ax[2].set_ylabel(r'$I(t)$')


    ax[0].set_title((r'$\theta_{{\mathrm{{old}}}} = {:.3f}, \theta_{{\mathrm{{new}}}} = {:.3f}, A_{{\mathrm{{perturb}}}} = {:.3f}$'
                     ).format(theta_old, theta_new, A_perturb))
    
    #--------------#
    #     Grid     #
    #--------------#
    for i in range(3):
        ax[i].grid(which='major')
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    # fig.savefig(filename_out)
    fig.show()

#------------------------------------------------------------------------------#
# Plot theta_new against min of I_seg4
def plot_theta_new_I_seg4(run_in):
    """
    Plots the minimum of the intensity of the perturbed segment against the
    returned phase theta_new.
    """
    import matplotlib.pyplot as plt
    from numpy import arange, mod, concatenate
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')
    
    #--------------------------------------------------------------------------#
    #                                Read Data                                 #
    #--------------------------------------------------------------------------#
    #-------------------------#
    #     Read Data: Data     #
    #-------------------------#
    # Read list of solutions
    solutions = run_in()

    # Empty arrays for data
    theta_new  = []
    I_seg4     = []

    # Cycle through solutions
    for i in range(len(solutions)):
    # Read solution
        sol = solutions[i]

        # Read theta_new
        theta_new_read = sol['theta_new']
        # Read seg4_x3
        I_seg4_read    = sol['seg4_x3']

        # Append data
        theta_new.append(theta_new_read)
        I_seg4.append(max(I_seg4_read))

    #--------------------------#
    #     Read: Parameters     #
    #--------------------------#
    # Read solution
    sol = solutions[0]

    # Periodicity
    k = sol['k']

    #--------------------------------------------------------------------------#
    #                                Plot Data                                 #
    #--------------------------------------------------------------------------#
    fignum = 'theta_new against I'
    plt.close(fignum)
    fig = plt.figure(num=fignum, figsize=[12, 8])
    ax  = plt.gca()
    
    #--------------#
    #     Plot     #
    #--------------#
    # Plot
    ax.plot(theta_new, I_seg4, color='C0', ls='solid')

    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # ax.set_xticks(arange(0.0, k+5, 5))
    # ax.set_xticks(arange(2.5, k+5, 5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    # ax.set_xlim(-0.01, k+0.1)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$\theta_{\mathrm{new}}$')
    ax.set_ylabel(r'$\mathrm{min} \left( I^{(4)} \right)$')

    # ax[0].set_title((r'$\theta_{{\mathrm{{old}}}} = {:.3f}, \theta_{{\mathrm{{new}}}} = {:.3f}, A_{{\mathrm{{perturb}}}} = {:.3f}$'
    #                  ).format(theta_old, theta_new, A_perturb))
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid(which='major')
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    # fig.savefig(filename_out)
    fig.show()