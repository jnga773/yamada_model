#==============================================================================#
#                    SELF-DEFINED FUNCTIONS FOR AUTO THINGS                    #
#==============================================================================#
# This Python module contains the following functions:
#
# - calc_stationary_points            : Calculates the three stationary
#                                       points of the Yamada model.
#
# - clean_directories                 : Deletes all of the AUTO generated
#                                       files that are yuck as.
#
# - save_move_data                    : Saves the data from run [RUN_IN] to
#                                       ./data/RUN_NAME_IN/.
#
# - bd_read                           : Loads the bifucation and solution data
#                                       in ./data_run_name_in/XXX.run_name_in.


#------------------------------------------------------------------------------#
def clean_directories():
    """
    Cleans up fort files, .o and .exe files in ./functions/, and the
    GLOBAL_VARIABLE .o file.
    """
    from os import remove, listdir

    # Cycle through files in current directory
    function_files = listdir('./')
    for file in function_files:
        # Remove c.XXX continuation scripts
        if file.startswith('c.'):
            remove('./{}'.format(file))
        
        # Remove Fortran .mod file
        if file.endswith('.mod'):
            remove('./{}'.format(file))
        
        # Remove fort.* data files
        if file.startswith('fort.'):
            remove('./{}'.format(file))


    # Remove things inside ./functions/ folder
    function_files = listdir('./functions/')
    for file in function_files:
        if file.endswith('.o') or file.endswith('.exe'):
            remove('./functions/{}'.format(file))

    print("Removed fort.*, *.o, *.exe, c.*, and *.mod files")
    
#------------------------------------------------------------------------------#
def save_move_data(run_in, run_name_in):
    """
    Saves the data from run [RUN_IN] to './data/RUN_NAME_IN/'.
    To make things easy, make sure that run_in and run_name_in are
    the same thing, except run_name_in will have '' around it :)
    """
    from os.path import isdir
    from os import makedirs
    import auto

    # Check for './data/' directory
    if not isdir('./data'):
        makedirs('./data')

    # Check for './data/run_name_in/' directory
    if not isdir('./data/{}'.format(run_name_in)):
        makedirs('./data/{}'.format(run_name_in))

    # Save the data
    auto.save(run_in, 'dat')

    # Move the data to the './data/run_name_in/' directory
    auto.move('dat', './data/{}'.format(run_name_in))

#------------------------------------------------------------------------------#
def bd_read(run_name_in):
    """
    Loads the bifurcation and solution data in ./data/run_name_in/XXX.run_name_in.
    """
    from os import remove, chdir
    from os.path import abspath
    import auto

    # Get current directory complete path
    main_dir = abspath('./')

    # Move into data directory
    chdir('./data/{}/'.format(run_name_in))

    # Load it
    bd_out = auto.loadbd('dat')

    # Change back to working directory
    chdir(main_dir)

    return bd_out
