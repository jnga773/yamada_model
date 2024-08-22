#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:54:58 2024

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('./figure_style.mplstyle')

plt.close('all')

# Figure filename
filename_out = "../../images/sect3/fig4_spectrum_labelled_high_low_drive.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def read_data():
    """
    Reads the PTC scan data from ./data/run09_PTC_scan/'
    """
    from auto import loadbd
    from PTC_plots import read_PTC_scan_data
    
    # Read data
    theta_old, theta_new, A_perturb = read_PTC_scan_data('run09_PTC_scan')
    
    # Read theta_perturb for d_vec components
    sol = 
    
    return theta_old, theta_new, A_perturb
    
    
#-----------------------------------------------------------------------------#
#                                  READ DATA                                  #
#-----------------------------------------------------------------------------#
theta_old, theta_new, A_perturb = read_data()


#-----------------------------------------------------------------------------#
#                                    PLOT                                     #
#-----------------------------------------------------------------------------#