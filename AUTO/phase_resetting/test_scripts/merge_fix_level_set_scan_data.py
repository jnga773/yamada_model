#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 16:05:30 2024

@author: jnga773
"""

import auto
import python_files.write_data_functions as data_funcs
import plotting_scripts.PTC_plots as plot_PTC

run_str_in = 'run09_level_set_scan'

# Get sorted data folders
data_dir = plot_PTC.sort_data_folders(run_str_in)

# Empty arrays for data
A_perturb     = []
theta_perturb = []
theta_new     = []

# Read data
for idx, run in enumerate(data_dir):
    # Read data
    bd = data_funcs.bd_read('{}/{}'.format(run_str_in, run))
    
    # # Split up bd runs
    # bd1 = bd[0]
    # bd2 = bd[1]

    # Merge runs
    bd_merge = auto.merge(bd)
    
    # Write and move data
    data_funcs.save_move_data(bd_merge, '{}/{}'.format(run_str_in, run))