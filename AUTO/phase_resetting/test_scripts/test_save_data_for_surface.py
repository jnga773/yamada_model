#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:44:45 2024

@author: jacob
"""

#--------------------------------------------------------------------------------------#
def sort_data_folders(run_str_in):
    from os import listdir
    from numpy import argsort
    
    # List all directories in ./data/run08_isochron_scan/
    dir_main = './data/{}/'.format(run_str_in)
    dir_list = listdir(dir_main)

    dir_list = [f for f in dir_list if not f.startswith('.')]
    
    # Cyycle through and get the integer number
    int_list = []
    for f in dir_list:
        int_read = f[4:]
        int_list.append(int(int_read))
        
    # Sort int_list and get indices
    idx_sorted = argsort(int_list)
    
    dir_list_sorted = [dir_list[i] for i in idx_sorted]
    
    return dir_list_sorted

#--------------------------------------------------------------------------------------#
def read_PTC_scan_bd_files(run_str_in):
    """
    Reads all of the b.dat files from the scan run and outputs a list of
    bd files.
    """
    from python_files.write_data_functions import bd_read

    # Get sorted data folders
    data_dir = sort_data_folders(run_str_in)
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    bd_out = []
    
    # Read data
    for idx, run in enumerate(data_dir):
        # Read bifurcation data
        bd = bd_read('{}/{}'.format(run_str_in, run))
        bd_out.append(bd)
        
    #----------------#
    #     Output     #
    #----------------#
    return bd_out

#--------------------------------------------------------------------------------------#
def check_runs_MX(bd_list_in):
    """
    Check which runs have MX in both directions, and which runs don't MX in
    either direction.
    """
    # Empty array for runs which MX
    MX_both = []
    not_MX_gt1 = []
    not_MX_lt1 = []
    
    for idx, bd in enumerate(bd_list_in):
        # Split bd runs
        bd_gt1 = bd[0]
        bd_lt1 = bd[1]
        
        # Check MX
        MX_gt1 = bd_gt1('MX')
        MX_lt1 = bd_lt1('MX')
        
        # Check if both branches MX
        if MX_gt1 and MX_lt1:
            # Both branches MX so save that label
            MX_both.append(idx)
        else:
            if not MX_gt1:
                # No MX in the gt1 branch
                not_MX_gt1.append(idx)
            if not MX_lt1:
                # No MX in the lt1 branch
                not_MX_lt1.append(idx)
                
    #----------------#
    #     Output     #
    #----------------#
    return MX_both, not_MX_gt1, not_MX_lt1            

#--------------------------------------------------------------------------------------#
def read_PTC_scan_data(bd_list_in):
    from numpy import flip

    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read a solution to get the unchanging parameters
    run_in   = bd_list_in[0]
    sol_read = run_in(1)

    # Perturbation direction angle
    theta_perturb = sol_read['theta_perturb']
    # System parameters
    gamma         = sol_read['gamma']
    A             = sol_read['A']
    B             = sol_read['B']
    a             = sol_read['a']
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Empty arrays for data
    theta_old_gt1 = []
    theta_old_lt1 = []
    theta_new_gt1 = []
    theta_new_lt1 = []
    A_perturb_out = []

    # Read data
    for idx, bd in enumerate(bd_list_in):
        # Get two directions
        bd_gt1 = bd[0]
        bd_lt1 = bd[1]

        #-----------------------#
        #     theta_old > 1     #
        #-----------------------#
        # Phase values
        old_gt1_read = bd_gt1['theta_old']
        new_gt1_read = bd_gt1['theta_new']

        #-----------------------#
        #     theta_old < 1     #
        #-----------------------#
        # Phase values
        old_lt1_read = flip(bd_lt1['theta_old'])
        new_lt1_read = flip(bd_lt1['theta_new'])

        #--------------------#
        #     Parameters     #
        #--------------------#
        sol = bd(1)
        A_perturb_read = sol['A_perturb']

        #----------------------#
        #      Append Data     #
        #----------------------#
        theta_old_gt1.append(old_gt1_read)
        theta_old_lt1.append(old_lt1_read)
        theta_new_gt1.append(new_gt1_read)
        theta_new_lt1.append(new_lt1_read)
        A_perturb_out.append(A_perturb_read)

    #----------------#
    #     Output     #
    #----------------#
    # Dictionary for parameters
    param_out = {'A': A, 'gamma': gamma, 'B': B, 'a': a,
                 'theta_perturb': theta_perturb}
    
    # Put into theta_old and theta_new arrays
    theta_old_out = {'gt1': theta_old_gt1, 'lt1': theta_old_lt1}
    theta_new_out = {'gt1': theta_new_gt1, 'lt1': theta_new_lt1}

    return param_out, theta_old_out, theta_new_out, A_perturb_out

#--------------------------------------------------------------------------------------#
def fix_theta_data(theta_old_in, theta_new_in):
    """
    Fixed the >1 and <1 data so that theta_old <= 1 and
    theta_new also sticks to the right places
    """
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Empty arrays
    theta_old_gt1 = []
    theta_old_lt1 = []
    theta_new_gt1 = []
    theta_new_lt1 = []

    # Fix up datafor i in range(len(theta_old)):
    for i in range(len(theta_old_in['gt1'])):
        # Theta_old values
        old_gt1 = theta_old_in['gt1'][i] - 1.0
        old_lt1 = theta_old_in['lt1'][i]
        # Theta_new values
        new_gt1 = theta_new_in['gt1'][i]
        new_lt1 = theta_new_in['lt1'][i]
        
        
        
        # print(f'old_gt1[0] = {old_gt1[0]:.3f} | old_gt1[-1] = {old_gt1[-1]:.3f}')
        # print(f'old_lt1[0] = {old_lt1[0]:.3f} | old_lt1[-1] = {old_lt1[-1]:.3f}')
        # print(f'new_gt1[0] = {new_gt1[0]:.3f} | new_gt1[-1] = {new_gt1[-1]:.3f}')
        # print(f'new_lt1[0] = {new_lt1[0]:.3f} | new_lt1[-1] = {new_lt1[-1]:.3f}')
        # print('\n')        
        
        # Check if theta_new goes below 0 at all
        if min(new_gt1) < 0.0:
            new_gt1 += 1.0
        if min(new_lt1) < 0.0:
            new_lt1 += 1.0
        
        # # Check if theta_new > 1 starts at 0.9 or higher. If so, shift
        # if min(new_gt1) > 0.85:
        #     new_gt1 += -1.0
        
        # Two data sets to plot
        theta_old_gt1.append(old_gt1)
        theta_old_lt1.append(old_lt1)
        
        theta_new_gt1.append(new_gt1)
        theta_new_lt1.append(new_lt1)

    #----------------#
    #     Output     #
    #----------------#
    # Put into theta_old and theta_new arrays
    theta_old_out = {'gt1': theta_old_gt1, 'lt1': theta_old_lt1}
    theta_new_out = {'gt1': theta_new_gt1, 'lt1': theta_new_lt1}
    
    return theta_old_out, theta_new_out

#--------------------------------------------------------------------------------------#
def pad_arrays(theta_old_in, theta_new_in, A_perturb_in):
    """
    Pads the theta_old and theta_new arrays to match the arrays with the max
    length. This makes things easier to plot using MATLAB's surf.
    """
    from numpy import array, pad, ones
    
    # Cycle through and calculate max length
    array_lengths = []
    for i in range(len(theta_old_in)):
        array_lengths.append(len(theta_old_in[i]))
    
    # Find maximum array length to pad to
    max_len = max(array_lengths)
    
    #------------------#
    #     Pad Data     #
    #------------------#
    # Pad array
    theta_old_out = []
    theta_new_out = []
    A_perturb_out = []
    
    for i in range(len(theta_old_in)):
        # Read arrays
        A_read = A_perturb_in[i]
        old_read = theta_old_in[i]
        new_read = theta_new_in[i]
        
        # Pad length
        pad_len = max_len - len(old_read)
        
        # Pad arrays        
        old_pad = pad(old_read, (pad_len, 0),
                      'constant', constant_values=(old_read[0], 0))
        
        new_pad = pad(new_read, (pad_len, 0),
                      'constant', constant_values=(new_read[0], 0))
        
        # times by ones()
        A_pad = A_read * ones(max_len)
        
        # Append to output
        theta_old_out.append(old_pad)
        theta_new_out.append(new_pad)
        A_perturb_out.append(A_pad)
    
    #----------------#
    #     Output     #
    #----------------#
    theta_old_out = array(theta_old_out).T
    theta_new_out = array(theta_new_out).T
    A_perturb_out = array(A_perturb_out).T
    
    return theta_old_out, theta_new_out, A_perturb_out

#------------------------------------------------------------------------------#
def save_PTC_scan(run_str_in):
    """
    Reads and saves the A_perturb, theta_old, and theta_new data
    from the scan run 'run09_PTC_scan' to a MATLAB .mat file.

    Input
    -------
    run_str_in : str
        String label identifier for the current run.
    """
    from scipy.io import savemat

    #-------------------#
    #     Read Data     #
    #-------------------#
    # Read all bifucation data files
    bd_list = read_PTC_scan_bd_files(run_str_in)
    
    # Initial read data
    param, theta_old_read, theta_new_read, A_perturb_read = read_PTC_scan_data(bd_list)
    
    # Get array of double MX
    MX_both, not_MX_gt1, not_MX_lt1 = check_runs_MX(bd_list)
    
    # Fix data
    theta_old_fix, theta_new_fix = fix_theta_data(theta_old_read, theta_new_read)
    
    # Pad data: theta_old >= 1
    old_gt1_pad, new_gt1_pad, A_perturb_gt1 = \
        pad_arrays(theta_old_fix['gt1'], theta_new_fix['gt1'], A_perturb_read)
        # Pad data: theta_old < 1
    old_lt1_pad, new_lt1_pad, A_perturb_lt1 = \
        pad_arrays(theta_old_fix['lt1'], theta_new_fix['lt1'], A_perturb_read)
    
    #-------------------#
    #     Save Data     #
    #-------------------#
    # Create dictionary for all saved stuff
    mat_out = {'A': param['A'], 'gamma': param['gamma'], 'B': param['B'], 'a': param['a'],
               'theta_perturb': param['theta_perturb'],
               'MX_both': MX_both, 'not_MX_gt1': not_MX_gt1, 'not_MX_lt1': not_MX_lt1,
               'theta_old_read': theta_old_read, 'theta_new_read': theta_new_read,
               'theta_old_fix': theta_old_fix, 'theta_new_fix': theta_new_fix,
               'theta_old_gt1': old_gt1_pad, 'theta_new_gt1': new_gt1_pad, 'A_perturb_gt1': A_perturb_gt1,
               'theta_old_lt1': old_lt1_pad, 'theta_new_lt1': new_lt1_pad, 'A_perturb_lt1': A_perturb_lt1}

    # Save data
    savemat('./data_mat/PTC_scan.mat', mat_out)
    
#=============================================================================#
if __name__ == '__main__':

    run_str = 'run09_PTC_scan_G'
    
    


    #-------------------#
    #     Save Data     #
    #-------------------#
    save_PTC_scan(run_str)
    