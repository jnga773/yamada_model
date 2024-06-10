#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 12:38:28 2024

@author: jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

def bd_read(run_name_in):
    """
    Loads the bifurcation and solution data in ./data/run_name_in/XXX.run_name_in.
    """
    from os import remove
    import auto

    # Copy the data to main directory
    auto.copy('./data/{}/'.format(run_name_in), 'dat')

    # Load it
    bd_out = auto.loadbd('dat')

    # Remove the copied data
    remove('./b.dat')
    remove('./s.dat')
    remove('./d.dat')

    return bd_out


#-------------------#
#     Read Data     #
#-------------------#
runs      = []
theta_old = []
theta_new = []
A_perturb = []

# Read data from bd files
for i in [2, 3, 9, 11]:
    # Run name
    run_name = 'run09_PTC_scan/sol_{}'.format(i)
    
    # Read data
    this_run = bd_read(run_name)
    # Append this run to runs array
    runs.append(this_run)
    
    # Read theta values
    theta_old_read = this_run['theta_old']
    theta_new_read = this_run['theta_new']
    # Append to arrays
    theta_old.append(theta_old_read)
    theta_new.append(theta_new_read)
    
    # Read A_perturb value
    A_perturb_read = this_run('UZ1')['A_perturb']
    # Append to array
    A_perturb.append(A_perturb_read)
    
# Fix data
for i in range(len(theta_old)):
    if min(theta_old[i]) > 0.8:
        theta_old[i] += -1.0
    if min(theta_new[i]) > 0.8:
        theta_new[i] += -1.0

#--------------#
#     Plot     #
#--------------#
fig = plt.figure(num='PTC Plot', figsize=[8, 8])
ax = plt.gca()

# Plot
for i in range(len(theta_old)):
    ax.plot(theta_old[i], theta_new[i], color='C{}'.format(i),
            label=r'$A_{{\mathrm{{perturb}}}} = {:.2f}$'.format(A_perturb[i]))
   
# Plot diagonal
ax.plot([0.0, 1.0], [0.0, 1.0], ls='dashed', color='k', alpha=0.75)

# Legend
ax.legend(loc='lower right')

# Ticks
ax.set_xticks(np.arange(0.0, 1.2, 0.2))
ax.set_xticks(np.arange(0.1, 1.2, 0.2), minor=True)

ax.set_yticks(np.arange(0.0, 1.2, 0.2))
ax.set_yticks(np.arange(0.1, 1.2, 0.2), minor=True)

# Limits
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(-0.01, 1.01)

# Grid
ax.grid()

# Axis titles
ax.set_xlabel(r'$\theta_{\mathrm{old}}$')
ax.set_ylabel(r'$\theta_{\mathrm{new}}$')

# Figure stuff
fig.tight_layout()
fig.savefig('./images/PTC_nonlinear.pdf')
fig.show()