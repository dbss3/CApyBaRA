###########################################################
import numpy as np
import sys
from read_in_routines import read_in_for_2DCARMA, read_file_to_dict
from check_files import *

###########################################################
# Aim here is to specify the location of a series of CARMA 
# output files from multiple restarts for 2DCARMA runs
# so that these can be read-in and then converted into 
# a single npz file that can be used for time-averaging etc

###########################################################
# Rough Outline:

# Already have function that reads in a CARMA routine and
# saves it, so run that multple times based on a 
# loc and path list

# Saving each in turn to npz has two advantages, can then
# still plot them individually if you wanted

# Then re-load the saved npz files and combine

# See if I can specify the times of each timestep

# Write a function that can then slice the combined struct
# based on a desired time window, i.e from t1 - t2

# Save this to a specific output npz, same structure as 
# a normal one though, so that it can be read in by your
# plotting routines.

###########################################################
def main():
    # read in file path for the file with the locs and path lists
    file_path = sys.argv[1]  # Get the file path from the command-line arguments

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

    # Need to make sure input file names and locations
    # have been provided
    # First the infile - i.e. the  2DCARMA output file
    infile_path = check_for_inpath(file_locs_and_names_dict,'infile')

    # Now the longitudes file, has to be read in seperately due
    # to no headers in the 2DCARMA output
    longitudes_path = check_for_inpath(file_locs_and_names_dict,'longitudes')

    # OK finally do the read in
    saved_dict_paths_list = read_in_for_2DCARMA(infile_path,longitudes_path,outfile_loc,run_name,\
                                            cloud_properties_file_path,cloud_materials_file_path)
