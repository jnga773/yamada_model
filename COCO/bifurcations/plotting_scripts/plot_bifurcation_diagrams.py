#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:00:12 2023

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add style sheet
plt.style.use('./style.py')

plt.close('all')

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
# Function to read data
# def read_data(filename_in):
#     # Read and output data from .txt files
#     from numpy import genfromtxt
    
#     # Append directory
#     filename_in = "../data/{}".format(filename_in)
    
#     # Read pump current data (A)
#     A_out = genfromtxt(filename_in, usecols=0, dtype=float)
    
#     # Read decay data
#     gamma_out = genfromtxt(filename_in, usecols=1, dtype=float)
    
#     return A_out, gamma_out

def read_data(bd_in, file_in):
    # Read data
    A_read     = bd_in['A_{}'.format(file_in)]
    gamma_read = bd_in['gamma_{}'.format(file_in)]
    
    # Reshape
    A_out     = A_read.reshape(A_read.shape[1])
    gamma_out = gamma_read.reshape(gamma_read.shape[1])
    
    return A_out, gamma_out
    
#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# # Hopf bifurcation data
# A_H, gamma_H = read_data('H.txt')

# # Neutral saddle=node
# A_SQ, gamma_SQ = read_data('SQ.txt')

# # Saddle-Node data
# A_S, gamma_S = read_data('S.txt')

# # Transcritical data
# A_T, gamma_T = read_data('T.txt')

# # Approximate homoclinic data
# A_L, gamma_L = read_data('L.txt')

# # Approximate homoclinic data
# A_D, gamma_D = read_data('D.txt')

# Read data matrix
from scipy.io import loadmat
bd = loadmat('../data/Bifurcation_Data.mat')

# Hopf bifurcation data
A_H, gamma_H = read_data(bd, 'H')

# Neutral saddle=node
A_SQ, gamma_SQ = read_data(bd, 'NSA')

# Saddle-Node data
A_S, gamma_S = read_data(bd, 'SN')

# Transcritical data
A_T, gamma_T = read_data(bd, 'T')

# Approximate homoclinic data
A_L, gamma_L = read_data(bd, 'homoclinic')

# Approximate homoclinic data
A_D, gamma_D = read_data(bd, 'double_limit')

#-----------------------------------------------------------------------------#
#                             PLOT (Full Picture)                             #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='Bifurcation Diagram', figsize=[12, 8])
ax = plt.gca()

# Plot data
ax.plot(A_S, gamma_S, color='C0', ls='solid', label='S')
ax.plot(A_T, gamma_T, color='C1', ls='dashed', label='T')
ax.plot(A_H, gamma_H, color='C2', ls='solid', label='H')
ax.plot(A_SQ, gamma_SQ, color='C3', ls='solid', label='SQ')
ax.plot(A_L, gamma_L, color='C4', ls='solid', label='L')
ax.plot(A_D, gamma_D, color='C5', ls='solid', label='D')

# Legend
ax.legend()

# Ticks - X axis
ax.set_xticks(np.arange(5.0, 12.0, 1.0))
ax.set_xticks(np.arange(5.2, 11.2, 0.2), minor=True)

# Ticks - Y Axis
ax.set_yticks(np.arange(0.0, 0.30, 0.05))
ax.set_yticks(np.arange(0.01, 0.25, 0.01), minor=True)

# Limits
ax.set_xlim(5.0, 11.0)
ax.set_ylim(0.0, 0.25)

# Grid
ax.grid()

# Labels
ax.set_xlabel(r'$A$')
ax.set_ylabel(r'$\gamma$')

# Title
ax.set_title(r'Yamada Model with $\left(B = 5.80, a = 1.80 \right)$')

# Figure stuff
fig.tight_layout(pad=0)
# fig.savefig("../images/Region_II_bifurcations_python.png", dpi=800)
# fig.savefig("../images/Region_II_bifurcations_python.pdf")
fig.show()

#-----------------------------------------------------------------------------#
#                           PLOT (Zoomed-In Picture)                          #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='Bifurcation Diagram (Zoomed)', figsize=[12, 8])
ax = plt.gca()

# Plot data
ax.plot(A_S, gamma_S, color='C0', ls='solid', label='S')
ax.plot(A_T, gamma_T, color='C1', ls='dashed', label='T')
ax.plot(A_H, gamma_H, color='C2', ls='solid', label='H')
ax.plot(A_SQ, gamma_SQ, color='C3', ls='solid', label='SQ')
ax.plot(A_L, gamma_L, color='C4', ls='solid', label='L')
ax.plot(A_D, gamma_D, color='C5', ls='solid', label='D')

# Straight line for A in Region 5
gamma_test = 0.035433
A_test = 7.3757
ax.axhline(y=gamma_test, color='k', alpha=0.5, label=rf'$\gamma = {gamma_test}$')
ax.axvline(x=A_test, color='k', alpha=0.5, label=rf'$A = {A_test}$')

# Legend
ax.legend()

# Ticks - X axis
ax.set_xticks(np.arange(6.6, 7.6, 0.2))
ax.set_xticks(np.arange(6.5, 7.45, 0.05), minor=True)

# Ticks - Y Axis
ax.set_yticks(np.arange(0.04, 0.11, 0.01))
ax.set_yticks(np.arange(0.042, 0.12, 0.002), minor=True)

# Limits
ax.set_xlim(6.5, 7.4)
ax.set_ylim(0.04, 0.10)

# Grid
ax.grid()

# Labels
ax.set_xlabel(r'$A$')
ax.set_ylabel(r'$\gamma$')

# Title
ax.set_title(r'Yamada Model with $\left(B = 5.80, a = 1.80 \right)$')

# Figure stuff
fig.tight_layout(pad=0)
# fig.savefig("../images/Region_II_bifurcations_zoomed_python.png", dpi=800)
# fig.savefig("../images/Region_II_bifurcations_zoomed_python.pdf")
fig.show()