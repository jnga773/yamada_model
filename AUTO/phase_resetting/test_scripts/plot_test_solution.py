#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:09:34 2024

@author: jnga773
"""

# Load extra functions
import auto
import python_files.write_data_functions as data_funcs
import plotting_scripts.PTC_plots as plot_PTC

# Load other stuff
import matplotlib.pyplot as plt
plt.close('all')

# Add thesis style sheet
plt.style.use('./plotting_scripts/figure_style.mplstyle')

# Load stable manifold
from scipy.io import loadmat
initial_PO = loadmat('./initial_PO.mat')

W_q = initial_PO['W_q_stable']

#-------------------#
#     Read Data     #
#-------------------#
run_name = 'run09_PTC_single'

bd = data_funcs.bd_read(run_name)

# Pick a solution
sol = bd(38)

# Initial orbit
x_PO = sol['seg1_x1']
y_PO = sol['seg1_x2']
z_PO = sol['seg1_x3']

# Perturb orbit
x_data = sol['seg4_x1']
y_data = sol['seg4_x2']
z_data = sol['seg4_x3']

#-------------------#
#     Plot Data     #
#-------------------#
fig = plt.figure(num='test', figsize=[12, 8])
ax = plt.axes(projection='3d')

# Plot initial orbit
ax.plot(x_PO, y_PO, z_PO, ls='solid', color='C2', lw=2.5, label=r'$\Gamma$')

ax.plot(x_data, y_data, z_data, ls='solid', lw=0.5, color='C0', alpha=0.75, label=r'$\Gamma_{p}$')

ax.plot(W_q[:, 0], W_q[:, 1], W_q[:, 2], ls='solid', label=r'$W^{s}(q)$')

# Legend
# ax.legend()

# Labels
# ax.set_xlabel(r'$G(t)$')
# ax.set_ylabel(r'$I(t)$')

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

# Figure stuff
# fig.tight_layout()
fig.show()