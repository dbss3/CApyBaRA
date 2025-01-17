import numpy as np
import os, sys
from read_in_routines import read_in_for_2DCARMA, read_file_to_dict

# Function that does everything
def main():
    """
    Main function to do everything
    """

    # Ensure exactly one argument is provided
    if len(sys.argv) != 3:
        print("Usage: python3 just_plot_it.py <file_path> <new_or_time_averaged_or_plotting_or_replotting=1,2,3,4>")
        sys.exit(1)

    # read in system argumnts
    file_path = sys.argv[1]  # Get the file path from the command-line arguments
    new_or_time_averaged_or_plotting_or_replotting = sys.argv[2]

    # Try reading in the file if it exists
    try:
        file_locs_and_names_dict = read_file_to_dict(file_path)
        print(f'File locations and names: {file_locs_and_names_dict}')
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    if new_or_time_averaged_or_plotting_or_replotting == 1:
        # This is a brand new attempt at plotting
        # So need to run the read in routine

        # Need to make sure input file names and locations
        # have been provided
        # First the infile - i.e. the  2DCARMA output file
        infile_path = check_for_path(file_locs_and_names_dict,'infile')

        # Now the longitudes file, has to be read in seperately due
        # to no headers in the 2DCARMA output
        longitudes_path = check_for_path(file_locs_and_names_dict,'longitudes')

        cloud_properties_file_path = check_for_path(file_locs_and_names_dict,'cloud_properties_file')

        cloud_materials_file_path = check_for_path(file_locs_and_names_dict,'cloud_materials_file')
        
        # OK finally do the read in
        if cloud_properties_file_path != 'use_default' and cloud_materials_file_path != 'use_default':
            saved_names_list = read_in_for_2DCARMA(infile_path,longitudes_path,infile_name,\
            cloud_properties_file_path=cloud_properties_file_path,
            cloud_materials_file_path=cloud_materials_file_path)
        elif cloud_properties_file_path == 'use_default' and cloud_materials_file_path != 'use_default':
            saved_names_list = read_in_for_2DCARMA(infile_path,longitudes_path,infile_name,\
            cloud_materials_file_path=cloud_materials_file_path)
        elif cloud_properties_file_path != 'use_default' and cloud_materials_file_path == 'use_default':
            saved_names_list = read_in_for_2DCARMA(infile_path,longitudes_path,infile_name,\
            cloud_properties_file_path=cloud_properties_file_path)
        else: # Using the default files
            saved_names_list = read_in_for_2DCARMA(infile_path,longitudes_path,infile_name)

    if new_or_time_averaged_or_plotting_or_replotting <= 2:
        # Have already read in once, so just need to
        # do time averaging
        pass

    if new_or_time_averaged_or_plotting_or_replotting <= 3:
        # Has already done the time averaging and 
        # just have to run the plotting scripts
        pass

    if new_or_time_averaged_or_plotting_or_replotting <= 4:
        # already done all the plotting before,
        # but probably want to change something
        # about the figure
        pass

    return

def check_for_path(file_locs_and_names_dict,file_str):
    if f'{file_str}_loc' in file_locs_and_names_dict:
        loc = file_locs_and_names_dict[f'{file_str}_loc']
    else:
        print(f'No {file_str}_loc for read in provided,')
        print('so assuming file is local to where the routine is being run (i.e. ./).')
        print('If this is not your desired effect,')
        print('please check your file_names_and_locs.txt file')
        loc = './'

    if f'{file_str}_name' in file_locs_and_names_dict:
        name = file_locs_and_names_dict[f'{file_str}_name']
    else:
        print(f'No {file_str}_name for read in provided,')

        # For criticial files, infile and longitudes
        # these MUST be provided so quit if not
        if file_str == 'infile' or file_str == 'longitudes':
            print('please check your file_names_and_locs.txt file')
            quit()
        else:
            print('However, this is a non-crucial file,')
            print('will default to the DEFAULT inputs.')
            path = 'use_default'

    # if crucial file, also check that the file exists
    if file_str == 'infile' or file_str == 'longitudes':

        path = f'{loc}/{name}.txt'

        if os.path.isfile(path):
            print(f"The file exists at: {path}")
            return path
        else:
            print(f"The file does not exist at: {path}")
            quit()
    
    return path

##########################################################
# Litle bit of python that means calling this script runs the function
if __name__ == "__main__":
    main()