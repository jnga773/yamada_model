#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 10:37:40 2025

@author: jnga773
"""

#------------------------------------------------------------------------------#
def plot_single_DTC(run_in):
    """
    Plots a single DTC

    Input
    -------
    run_in : AUTO generated run file
        The continuation rum.
    """
    import matplotlib.pyplot as plt
    from numpy import arange, cos, sin, pi
    import matplotlib.patches as patches
    
    # Add thesis style sheet
    plt.style.use('./plotting_scripts/figure_style.mplstyle')

    # # Filename for figure
    # filename_out = "./images/PTC_single.pdf"

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Cooridnates of the isochron
    theta_perturb = run_in['theta_perturb']
    theta_new     = run_in['theta_new']
    
    # Sort for theta_perturb < 0 and theta_perturb > 0
    mask_lt0 = theta_perturb < 0
    mask_gt0 = theta_perturb >= 0
    
    theta_perturb_lt0 = theta_perturb[mask_lt0] + 1.0
    theta_perturb_gt0 = theta_perturb[mask_gt0]
    theta_new_lt0     = theta_new[mask_lt0]
    theta_new_gt0     = theta_new[mask_gt0]
    

    # Solution to read
    sol_read = run_in(1)

    # Perturbation amplitude
    A_perturb     = sol_read['A_perturb']

    #-------------------------------------------------------------------------#
    #                                   Plot                                  #
    #-------------------------------------------------------------------------#
    fig = plt.figure(num='Single DTC Test', figsize=[8, 8])
    fig.clear()
    ax = plt.gca()

    #--------------#
    #     Plot     #
    #--------------#
    # Plot fundamental domain
    rect = patches.Rectangle((-1, 0), 2, 1, edgecolor='none', facecolor='C2', alpha=0.25)
    ax.add_patch(rect)
    
    # Plot DTC (theta_perturb >= 0)
    ax.plot(theta_perturb_gt0, theta_new_gt0, color='C0', ls='solid',
            label=r'$\varphi_{\bm{d}} \geq 0$')
    # Plot DTC (theta_perturb < 0)
    ax.plot(theta_perturb_lt0, theta_new_lt0, color='C1', ls='dashed',
            label=r'$\varphi_{\bm{d}} < 0$')
    
    # Legend
    ax.legend(loc='upper right')
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # ax.set_xticks(arange(0.0, 1.2, 0.2))
    # ax.set_xticks(arange(0.1, 1.2, 0.2), minor=True)
    
    # ax.set_yticks(arange(0.0, 3.5, 0.5))
    # ax.set_yticks(arange(0.25, 3.5, 0.5), minor=True)
    
    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax.set_xlim(0, 1)
    # ax.set_ylim(-0.1, 1.1)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax.set_xlabel(r'$\varphi_{\bm{d}}$')
    ax.set_ylabel(r'$\vartheta_{\mathrm{n}}$')
    ax.set_title((r'Single DTC (TEST) with $A = {:.3f}$'
                  ).format(A_perturb))
    
    #--------------#
    #     Grid     #
    #--------------#
    ax.grid()
    
    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    fig.tight_layout()
    # fig.savefig(filename_out)
    fig.show()