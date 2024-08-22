#==============================================================================#
#                  FUNCTIONS USED IN THE PTC SCAN CALCULATION                  #
#==============================================================================#
# Bunch of functions to save the PTC scan data

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
    from python_files.data_functions import bd_read

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
        
        # Check if theta_new goes below 0 at all
        if min(new_gt1) < 0.3:
            new_gt1 += 1.0
        if min(new_lt1) < 0.3:
            new_lt1 += 1.0
        
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

#------------------------------------------------------------------------------#
def sort_data_holes(theta_old_in, theta_new_in, bd_list_in, A_perturb_in, which):
    """
    Sorts the data into three sections: before the hole, the hole, and after
    the hole.
    which = {'hole', 'before', 'after'}
    """
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Find MX points in data
    MX_both, not_MX_gt1, not_MX_lt1 = check_runs_MX(bd_list_in)
    
    # Get 'fixed' data
    old_lt1 = theta_old_in['lt1']
    old_gt1 = theta_old_in['gt1']
    new_lt1 = theta_new_in['lt1']
    new_gt1 = theta_new_in['gt1']
    
    #-------------------------#
    #     Sort Data: Hole     #
    #-------------------------#
    if which == 'hole':
        # theta_old data in hole
        hole_old_gt1 = [old_gt1[i] for i in MX_both]
        hole_old_lt1 = [old_lt1[i] for i in MX_both]
    
        # theta_new data in hole
        hole_new_gt1 = [new_gt1[i] for i in MX_both]
        hole_new_lt1 = [new_lt1[i] for i in MX_both]
    
        # A_perturb data in hole
        hole_A = [A_perturb_in[i] for i in MX_both]
        
    #---------------------------#
    #     Sort: Before Hole     #
    #---------------------------#
    if which == 'before':
        # Empty arrays
        before_hole_old    = []
        before_hole_new    = []
        before_hole_labels = []
    
        for i in range(len(A_perturb_in)):
            # Before the hole
            if i < MX_both[0]:
                # Save indices
                before_hole_labels.append(i)
                
                if i in not_MX_gt1 and i in not_MX_lt1:
                    # print(str(i).zfill(3), 'full runs in both gt1 and lt1')
                    # Full runs in both directions, so append just the 'gt1' run
                    before_hole_old.append(old_lt1[i])
                    before_hole_new.append(new_lt1[i])
    
                else:
                    # Check if full run in 'gt1'
                    if i in not_MX_gt1:
                        # print(str(i).zfill(3), 'full run in gt1')
                        before_hole_old.append(old_gt1[i])
                        before_hole_new.append(new_gt1[i])
                    elif i in not_MX_lt1:
                        # print(str(i).zfill(3), 'full run in lt1')
                        before_hole_old.append(old_lt1[i])
                        before_hole_new.append(new_lt1[i])
        
        # A_perturb values
        before_hole_A = [A_perturb_in[i] for i in before_hole_labels]
        
    #---------------------------#
    #     Sort: After Hole     #
    #---------------------------#
    if which == 'after':
        # Empty arrays
        after_hole_old    = []
        after_hole_new    = []
        after_hole_labels = []
    
        for i in range(len(A_perturb_in)):
            # Before the hole
            if i > MX_both[-1]:
                # Save indices
                after_hole_labels.append(i)
                
                if i in not_MX_gt1 and i in not_MX_lt1:
                    # print(str(i).zfill(3), 'full runs in both gt1 and lt1')
                    # Full runs in both directions, so append just the 'gt1' run
                    after_hole_old.append(old_lt1[i])
                    after_hole_new.append(new_lt1[i])
    
                else:
                    # Check if full run in 'gt1'
                    if i in not_MX_gt1:
                        # print(str(i).zfill(3), 'full run in gt1')
                        after_hole_old.append(old_gt1[i])
                        after_hole_new.append(new_gt1[i])
                    elif i in not_MX_lt1:
                        # print(str(i).zfill(3), 'full run in lt1')
                        after_hole_old.append(old_lt1[i])
                        after_hole_new.append(new_lt1[i])
                        
        # A_perturb values
        after_hole_A = [A_perturb_in[i] for i in after_hole_labels] 
        
    #----------------#
    #     Output     #
    #----------------#
    if which == 'hole':
        data_hole_gt1 = {'theta_old': hole_old_gt1,
                         'theta_new': hole_new_gt1,
                         'A_perturb': hole_A}
        data_hole_lt1 = {'theta_old': hole_old_lt1,
                         'theta_new': hole_new_lt1,
                         'A_perturb': hole_A}
        
        return data_hole_gt1, data_hole_lt1
    
    elif which == 'before':
        data_before_hole = {'theta_old': before_hole_old,
                            'theta_new': before_hole_new,
                            'A_perturb': before_hole_A}
        
        return data_before_hole
        
    elif which == 'after':
        data_after_hole = {'theta_old': after_hole_old,
                           'theta_new': after_hole_new,
                           'A_perturb': after_hole_A} 
        
        return data_after_hole
   
#--------------------------------------------------------------------------------------# 
def closest(list_in, K_in):
    """
    Find closest value in list to K
    """
    from numpy import asarray, abs
    # Turn into array
    list_in = asarray(list_in)
    # Find point in list closest to K
    diff = abs(list_in - K_in)
    idx_out = diff.argmin()
    
    return idx_out

def fill_surface_gaps(data_hole_gt1_in, data_hole_lt1_in, data_before_hole_in, data_after_hole_in):
    """
    Merges the edges of the hole data to the not-hole data
    """
    #-------------------#
    #     Read Data     #
    #-------------------#
    # hole data (theta_old > 1)
    hole_old_gt1 = data_hole_gt1_in['theta_old']
    hole_new_gt1 = data_hole_gt1_in['theta_new']
    hole_A_gt1   = data_hole_gt1_in['A_perturb']
    
    # hole data (theta_old < 1)
    hole_old_lt1 = data_hole_lt1_in['theta_old']
    hole_new_lt1 = data_hole_lt1_in['theta_new']
    hole_A_lt1   = data_hole_lt1_in['A_perturb']
    
    # Before hole data
    before_old = data_before_hole_in['theta_old']
    before_new = data_before_hole_in['theta_new']
    before_A   = data_before_hole_in['A_perturb']
    
    # After hole data
    after_old = data_after_hole_in['theta_old']
    after_new = data_after_hole_in['theta_new']
    after_A   = data_after_hole_in['A_perturb']
    
    #-----------------------#
    #     Create Copies     #
    #-----------------------#
    # hole data (theta_old > 1)
    hole_old_gt1_merge = hole_old_gt1.copy()
    hole_new_gt1_merge = hole_new_gt1.copy()
    hole_A_gt1_merge   = hole_A_gt1.copy()
    
    # hole data (theta_old < 1)
    hole_old_lt1_merge = hole_old_lt1.copy()
    hole_new_lt1_merge = hole_new_lt1.copy()
    hole_A_lt1_merge   = hole_A_lt1.copy()
    
    #---------------------------------------#
    #    Merge Data: Before to hole_gt1     #
    #---------------------------------------#
    # Find point in (not_old) closest to the end point of hole_old
    idx_min = closest(before_old[-1], (hole_old_gt1[0])[-1])
    
    # Append truncated data to 'hole' data
    hole_old_gt1_merge.insert(0, (before_old[-1])[:idx_min])
    hole_new_gt1_merge.insert(0, (before_new[-1])[:idx_min])
    hole_A_gt1_merge.insert(0, before_A[-1])
    
    #---------------------------------------#
    #    Merge Data: Before to hole_lt1     #
    #---------------------------------------#
    # Find point in (not_old) closest to the end point of hole_old
    # idx_min = closest(before_old[-1], (hole_old_lt1[0])[-1])
    
    # Append truncated data to 'hole' data
    # hole_old_lt1_merge.insert(0, (before_old[-1])[idx_min:])
    # hole_new_lt1_merge.insert(0, (before_new[-1])[idx_min:]-1.0)
    hole_old_lt1_merge.insert(0, hole_old_lt1[0])
    hole_new_lt1_merge.insert(0, hole_new_lt1[0])
    hole_A_lt1_merge.insert(0, before_A[-1])
    
    #---------------------------------------#
    #     Merge Data: hole_gt1 to after     #
    #---------------------------------------#
    # Find point in (not_old) closest to the end point of hole_old
    idx_min = closest(after_old[0], (hole_old_gt1[-1])[-1])
    
    # Append truncated data to 'hole' data
    hole_old_gt1_merge.append((after_old[0])[:idx_min])
    hole_new_gt1_merge.append((after_new[0])[:idx_min])
    hole_A_gt1_merge.append(after_A[0])
    
    #---------------------------------------#
    #     Merge Data: hole_lt1 to after     #
    #---------------------------------------#
    # Find point in (not_old) closest to the end point of hole_old
    idx_min = closest(after_old[0], (hole_old_lt1[-1])[0])
    
    # Append truncated data to 'hole' data
    hole_old_lt1_merge.append((after_old[0])[idx_min:])
    hole_new_lt1_merge.append((after_new[0])[idx_min:])
    hole_A_lt1_merge.append(after_A[0])
    
    #----------------#
    #     Output     #
    #----------------#
    # theta_old > 1
    data_hole_gt1_out = {'theta_old': hole_old_gt1_merge,
                         'theta_new': hole_new_gt1_merge,
                         'A_perturb': hole_A_gt1_merge}
    # theta_old < 1
    data_hole_lt1_out = {'theta_old': hole_old_lt1_merge,
                         'theta_new': hole_new_lt1_merge,
                         'A_perturb': hole_A_lt1_merge}
    
    return data_hole_gt1_out, data_hole_lt1_out 

#--------------------------------------------------------------------------------------#
def pad_arrays(data_in):
    """
    Pads the theta_old and theta_new arrays to match the arrays with the max
    length. This makes things easier to plot using MATLAB's surf.
    """
    from numpy import array, pad, ones
    
    #---------------#
    #     Input     #
    #---------------#
    theta_old_in = data_in['theta_old']
    theta_new_in = data_in['theta_new']
    A_perturb_in = data_in['A_perturb']
    
    #-------------------------#
    #     Find Max Length     #
    #-------------------------#    
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
    
    # Turn into a data structure / dictionary
    data_out = {'theta_old': theta_old_out,
                'theta_new': theta_new_out,
                'A_perturb': A_perturb_out}
    
    return data_out

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
    
    #-------------------#
    #     Sort Data     #
    #-------------------#
    # Fix data
    theta_old_fix, theta_new_fix = fix_theta_data(theta_old_read, theta_new_read)
    
    # Sort data inside hole
    data_hole_gt1, data_hole_lt1 = sort_data_holes(theta_old_fix, theta_new_fix, bd_list, A_perturb_read, 'hole')

    # Sort data before hole
    data_before_hole = sort_data_holes(theta_old_fix, theta_new_fix, bd_list, A_perturb_read, 'before')

    # Sort data after hole
    data_after_hole = sort_data_holes(theta_old_fix, theta_new_fix, bd_list, A_perturb_read, 'after')
    
    #---------------------------------------#
    #     Fill In Gaps Between Segments     #
    #---------------------------------------#
    # Merge hole_gt1 to before
    data_hole_gt1_merge, data_hole_lt1_merge = fill_surface_gaps(data_hole_gt1, data_hole_lt1, data_before_hole, data_after_hole)

    #------------------#
    #     Pad Data     #
    #------------------#
    # Inside hole: theta_old > 1
    # data_hole_gt1_pad = pad_arrays(data_hole_gt1)
    data_hole_gt1_pad = pad_arrays(data_hole_gt1_merge)

    # Inside hole: theta_old < 1
    # data_hole_lt1_pad = pad_arrays(data_hole_lt1)
    data_hole_lt1_pad = pad_arrays(data_hole_lt1_merge)

    # Before hole
    data_before_hole_pad = pad_arrays(data_before_hole)

    # After hole
    data_after_hole_pad = pad_arrays(data_after_hole)
    
    #-------------------#
    #     Save Data     #
    #-------------------#
    # Create dictionary for all saved stuff
    mat_out = {'A': param['A'], 'gamma': param['gamma'], 'B': param['B'], 'a': param['a'],
               'theta_perturb': param['theta_perturb'], 'A_perturb': A_perturb_read}
    mat_out['data_hole_gt1']    = data_hole_gt1_pad
    mat_out['data_hole_lt1']    = data_hole_lt1_pad
    mat_out['data_before_hole'] = data_before_hole_pad
    mat_out['data_after_hole']  = data_after_hole_pad

    # Save data
    savemat('./data_mat/PTC_scan.mat', mat_out)