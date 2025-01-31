import numpy as np
import os, sys
import glob
from read_in_routines import read_in_for_2DCARMA, read_file_to_dict, read_cloud_properties
from time_averaging import do_time_averaging
from plotting_routines import plotter
from check_files import *

########################################################
# Function that does everything
def main():
    """
    Main function to do EVERYTHING
    """

    # Ensure exactly one argument is provided
    if len(sys.argv) != 5:
        print("Usage: python3 just_plot_it.py <file_path> <new_or_time_averaged_or_plotting_or_replotting=1,2,3,4>  <which_plot> <desire_value>")
        sys.exit(1)

    # read in system argumnts
    file_path = sys.argv[1]  # Get the file path from the command-line arguments
    new_or_time_averaged_or_plotting_or_replotting = int(sys.argv[2])
    which_plot = sys.argv[3]
    desire_value = float(sys.argv[4]) # MUST pass dummy value if the which_plot doesn't need one
    #TODO: will add so that more values can be passed

    # Try reading in the file if it exists
    try:
        file_locs_and_names_dict = read_file_to_dict(file_path)
        print(f'File locations and names: {file_locs_and_names_dict}\n')
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    #See if the output file is established
    outfile_loc = check_for_outpath(file_locs_and_names_dict,'outfile')

    #Allow the user to establish a different run_name, otherwise just use the infile_name
    run_name = check_for_run_name(file_locs_and_names_dict,'run_name')

    cloud_properties_file_path = check_for_inpath(file_locs_and_names_dict,'cloud_properties_file')

    cloud_materials_file_path = check_for_inpath(file_locs_and_names_dict,'cloud_materials_file')
        
    if new_or_time_averaged_or_plotting_or_replotting == 1:
        # This is a brand new attempt at plotting
        # So need to run the read in routine

        print('Option 1 start')

        # Need to make sure input file names and locations
        # have been provided
        # First the infile - i.e. the  2DCARMA output file
        infile_path = check_for_inpath(file_locs_and_names_dict,'infile')

        # Now the longitudes file, has to be read in seperately due
        # to no headers in the 2DCARMA output
        longitudes_path = check_for_inpath(file_locs_and_names_dict,'longitudes')

        # OK finally do the read in
        if type(infile_path) == str:
            saved_dict_paths_list = read_in_for_2DCARMA(infile_path,longitudes_path,outfile_loc,run_name,\
                                                        cloud_properties_file_path,cloud_materials_file_path)
        else: # assume list of infile_paths
            infile_path_list = infile_path
            saved_dict_paths_list = []
            for infile_path in infile_path_list:
                saved_dict_paths = read_in_for_2DCARMA(infile_path,longitudes_path,outfile_loc,run_name,\
                                                        cloud_properties_file_path,cloud_materials_file_path)
                saved_dict_paths_list.append(saved_dict_paths)

                # COMPILATOR: NEED FUNCTION HERE TO TURN THESE MULTIPLE 
                # SAVED DICTS INTO A SINGLE ONE, SO THE REST OF THE SCRIPT
                # IS UNCHANGED

        print('Option 1 end\n')

    if new_or_time_averaged_or_plotting_or_replotting <= 2:
        # Have just read in so saved_dict_paths_list is defined

        print('Option 2 start')

        if new_or_time_averaged_or_plotting_or_replotting == 1:
            # do time averaging, this loads the neccessary dicts and averages them
            pass
        else: # Need to go get the files that MUST exist already
            # Have already read in once, so just need to get the paths from the outfile_loc 
            # glob it up time Mr globby
            saved_dict_paths_list = glob.glob(f'{outfile_loc}/temporal_*.npz')
            #sort it as we want in alphabetical oreder (global_cloud_prop, mmr, svp)
            saved_dict_paths_list.sort()
        
        #load in the dictionaries
        loaded_dicts_list = []
        for saved_dict_path in saved_dict_paths_list:
            loaded_dict = np.load(saved_dict_path)
            loaded_dicts_list.append(loaded_dict)
        
        saved_time_averaged_dict_paths_list = do_time_averaging(loaded_dicts_list,outfile_loc,run_name)

        print('Option 2 end\n')

    if new_or_time_averaged_or_plotting_or_replotting <= 3:
        # Has already done the time averaging and saved_dict_paths are defined
        # just have to run the plotting scripts
        print('Option 3 start')

        if new_or_time_averaged_or_plotting_or_replotting <= 2:
            pass
        else: #everything up to time averaging done in the past
            #glob me up scotty some time averaged files
            saved_time_averaged_dict_paths_list = glob.glob(f'{outfile_loc}/time_averaged_*.npz')
            saved_time_averaged_dict_paths_list.sort()

        # for now I'm going to make it always plot all groups seperately
        # TODO: also add so that when it does the loop it makes the combined figure
        cloud_properties_dict = read_cloud_properties(cloud_properties_file_path)
        group_name_list = cloud_properties_dict['group_name']

        # For now only cloud property plots are supported
        # TODO: Add an option here based on the which_plot
        loaded_dict = np.load(saved_time_averaged_dict_paths_list[0])

        groups_with_errors = []
        for group_name in group_name_list:
            try:
                plotter(loaded_dict,group_name,which_plot,desire_value,outfile_loc,run_name,cloud_properties_dict)
            except:
                print(f'Error with plotting group: {group_name}')
                groups_with_errors.append(group_name)

        print('Finished running script. Errors with plotting:')
        for e in groups_with_errors:
            print(e)

        print('Option 3 end\n')

    if new_or_time_averaged_or_plotting_or_replotting <= 4:
        # already done all the plotting before,
        # but probably want to change something
        # about the figure
        #TODO: replotter function
        print('Option 4 start')

        print('Option 4 end\n')

    return

##########################################################
# Litle bit of python that means calling this script runs the function
if __name__ == "__main__":
    main()