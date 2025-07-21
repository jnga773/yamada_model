#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:28:47 2025

@author: jnga773
"""

from scipy.io import loadmat, savemat
import continuation_scripts.phase_reset.save_PTC_scan_data as shiet


# Load data
mat_out = loadmat('./data_mat/PTC_scan.mat')


data_hole_gt1 = mat_out['data_hole_gt1']
data_hole_lt1 = mat_out['data_hole_lt1']

data_before_hole = mat_out['data_before_hole']
data_after_hole = mat_out['data_after_hole']

theta_old_gt1 = mat_out['theta_old_gt1']
theta_old_lt1 = mat_out['theta_old_lt1']
theta_new_gt1 = mat_out['theta_new_gt1']
theta_new_lt1 = mat_out['theta_new_lt1']

theta_old_read = {'gt1': theta_old_gt1, 'lt1': theta_old_lt1}
theta_new_read = {'gt1': theta_new_gt1, 'lt1': theta_new_lt1}
A_perturb_read = mat_out['A_perturb']

#-------------------#
#     Sort Data     #
#-------------------#
# Fix data
theta_old_fix, theta_new_fix = \
    shiet.fix_theta_data(theta_old_read, theta_new_read)

# Sort data inside hole
data_hole_gt1, data_hole_lt1 = \
    shiet.sort_data_holes(theta_old_fix, theta_new_fix, bd_list, A_perturb_read, 'hole')

# Sort data before hole
data_before_hole = sort_data_holes(theta_old_fix, theta_new_fix, bd_list, A_perturb_read, 'before')

# Sort data after hole
data_after_hole = sort_data_holes(theta_old_fix, theta_new_fix, bd_list, A_perturb_read, 'after')

#---------------------------------------#
#     Fill In Gaps Between Segments     #
#---------------------------------------#
# Merge hole_gt1 to before
data_hole_gt1_merge, data_hole_lt1_merge = \
    shiet.fill_surface_gaps(data_hole_gt1, data_hole_lt1, data_before_hole, data_after_hole)