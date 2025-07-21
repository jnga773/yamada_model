#==============================================================================#
#                FUNCTIONS USED IN THE ISOCHRON SCAN CALCULATION               #
#==============================================================================#
# Bunch of functions to save the isochron scan data

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
def read_isochron_scan_bd_files(run_str_in):
    """
    Reads all of the b.dat files from the scan run and outputs a list of
    bd files.
    """
    from continuation_scripts.data_functions import bd_read

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
def read_isochron_scan_data(bd_list_in):
    #-------------------------#
    #     Read Parameters     #
    #-------------------------#
    # Read a solution to get the unchanging parameters
    run_in   = bd_list_in[0]
    sol_read = run_in(1)

    # Phase of periodic orbit
    theta_old = sol_read['theta_old']
    theta_new = sol_read['theta_new']
    # System parameters
    gamma     = sol_read['gamma']
    A         = sol_read['A']
    B         = sol_read['B']
    a         = sol_read['a']
    
    #-------------------#
    #     Read Data     #
    #-------------------#
    # Empty arrays for data
    iso1_data = []
    iso2_data = []
    iso3_data = []

    # Read data
    for idx, bd in enumerate(bd_list_in):
        #-------------------------#
        #     Isochron Values     #
        #-------------------------#
        # iso1
        iso1_read = bd['iso1']
        # iso2
        iso2_read = bd['iso2']
        # iso3
        iso3_read = bd['iso3']

        #----------------------#
        #      Append Data     #
        #----------------------#
        iso1_data.append(iso1_read)
        iso2_data.append(iso2_read)
        iso3_data.append(iso3_read)

    #----------------#
    #     Output     #
    #----------------#
    # Dictionary for parameters
    param_out = {'A': A, 'gamma': gamma, 'B': B, 'a': a,
                 'theta_old': theta_old, 'theta_new': theta_new}

    return param_out, iso1_data, iso2_data, iso3_data

#--------------------------------------------------------------------------------------#
def save_isochron_scan_data(run_str_in):
    #-------------------#
    #     Read Data     #
    #-------------------#
    # List of bd files
    bd_list = read_isochron_scan_bd_files(run_str_in)

    # Read data
    param_out, iso1_data, iso2_data, iso3_data = read_isochron_scan_data(bd_list)

    #-------------------#
    #     Save Data     #
    #-------------------#
    # Save data
    from scipy.io import savemat

    # Create dictionary for all saved stuff
    mat_out = {'parameters': param_out, 'iso1_data': iso1_data, 'iso2_data': iso2_data, 'iso3_data': iso3_data}

    # Save data
    savemat('./data_mat/isochron_scan.mat', mat_out)