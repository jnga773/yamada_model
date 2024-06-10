#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:10:59 2024

@author: jnga773
"""

#------------------------------------------------------------------------------#
# Plot the phase-space solution from the Hopf-to-PO continuation
def plot_hopf_to_PO_solution(sol_PO_in):
    """
    Plots the phase-space solution of the Hopf-to-PO continuation run.
    
    Input
    -------
    sol_PO_in : AUTO-continuation solution
        AUTO generate solution for the periodic orbit (like run('UZ')).
        Probably run03_hopf_to_PO('UZ1').
    """
    import matplotlib.pyplot as plt
    from numpy import arange
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')
    
    # Figure name for file
    filename_out = './images/hopf_to_PO_solution.pdf'
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    #------------------------#
    #     Periodic Orbit     #
    #------------------------#
    # Get state space data
    x1 = sol_PO_in['x1']
    x2 = sol_PO_in['x2']
    x3 = sol_PO_in['x3']
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Hopf to PO solution', figsize=[6, 6])
    ax = plt.axes(projection='3d')
    
    #--------------#
    #     Plot     #
    #--------------#
    # Plot periodic orbit
    ax.plot(x1, x2, x3, color='C0', ls='solid', label=r'$\Gamma$')
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax.set_xticks(arange(0.0, 5.5, 0.5))
    ax.set_xticks(arange(0.25, 5.5, 0.5), minor=True)
    
    ax.set_yticks(arange(0.0, 4.0, 0.5))
    ax.set_yticks(arange(0.25, 4.0, 0.5), minor=True)

    ax.set_zticks(arange(0.0, 22, 4.0))
    ax.set_zticks(arange(2.0, 22, 4.0), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0.0, 5.0)
    ax.set_ylim(0.0, 3.5)
    ax.set_zlim(0.0, 18)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$G(t)$')
    ax.set_ylabel(r'$Q(t)$')
    ax.set_zlabel(r'$I(t)$')

    ax.set_title(r'Initial Solution from Hopf Bifurcation')
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid(visible=True, which='major')
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    fig.savefig(filename_out)
    fig.show()   

#------------------------------------------------------------------------------#
# Plot initial phase space template
def plot_base_periodic_orbit(ax_in):
    """
    Plots the initial phase space portrait of the periodic orbit from the saved
    data in ./data/.

    Input
    -------
    ax_in : matplotlib.pyplot ax structure
        Input axis to plot the solutions.

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
    data = loadmat('./data/initial_PO.mat')

    # Plot things
    xbp  = data['xbp_PO']
    x0   = data['x0']
    xpos = data['xpos']
    xneg = data['xneg']

    #--------------#
    #     Plot     #
    #--------------#
    # Plot periodic orbit
    ax_in.plot(xbp[:, 0], xbp[:, 1], xbp[:, 2], color='C2', ls='solid', label=r'$\Gamma$')
    
    # Plot x0
    ax_in.plot(x0[0], x0[1], x0[2], color='r', ls='none', label=r'$o$',
              marker='o', markeredgecolor='r', markerfacecolor='r',
              markersize=12)
    
    # Plot xpos
    ax_in.plot(xpos[0], xpos[1], xpos[2], color='b', ls='none', label=r'$q$',
              marker='*', markeredgecolor='b', markerfacecolor='b',
              markersize=12)
    
    # Plot xneg
    ax_in.plot(xneg[0], xneg[1], xneg[2], color='r', ls='none', label=r'$p$',
              marker='*', markeredgecolor='r', markerfacecolor='r',
              markersize=12)
    
    #----------------#
    #     Output     #
    #----------------#
    return ax_in

#------------------------------------------------------------------------------#
# Plotting function for the shifted periodic orbit solution
def plot_initial_PO():
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
    filename_out = './images/initial_PO_solition.pdf'
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    # Load data from matrix
    data = loadmat('./data/initial_PO.mat')

    # Plot things
    xbp = data['xbp_PO']
    x0  = data['x0']
    
    #-------------------------------------------------------------------------#
    #                                Read Data                                #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Shifted Periodic Orbit', figsize=[6, 6])
    ax = plt.axes(projection='3d')
    
    #--------------#
    #     Plot     #
    #--------------#
    ax = plot_base_periodic_orbit(ax)
    
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
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, 8.0)
    ax.set_zlim(0.0, 20)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$G(t)$')
    ax.set_ylabel(r'$Q(t)$')
    ax.set_zlabel(r'$I(t)$')

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