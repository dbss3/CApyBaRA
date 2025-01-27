import os

########################################################
# Functions that check if given input files exist and/or
# are specified correctly

def check_for_inpath(file_locs_and_names_dict,file_str):
    if f'{file_str}_loc' in file_locs_and_names_dict:
        loc = file_locs_and_names_dict[f'{file_str}_loc']
    else:
        print(f'No {file_str}_loc for read in provided,')
        print('so assuming file is local to where the routine is being run (i.e. ./).')
        print('If this is not your desired effect,')
        print('please check your file_names_and_locs.txt file\n')
        loc = './'

    if f'{file_str}_name' in file_locs_and_names_dict:
        name = file_locs_and_names_dict[f'{file_str}_name']

    else:
        print(f'No {file_str}_name for read in provided,')

        # For criticial files, infile and longitudes
        # these MUST be provided so quit if not
        if file_str == 'infile' or file_str == 'longitudes':
            print('please check your file_names_and_locs.txt file\n')
            quit()
        elif file_str == 'cloud_properties_file':
            print('however, this is a non-crucial file,')
            print('will default to the DEFAULT inputs.')
            path = './group_names_and_properties_DEFAULT.txt'
            print(f'{path}\n')
        elif file_str == 'cloud_materials_file':
            print('however, this is a non-crucial file,')
            print('will default to the DEFAULT inputs.')
            path = './cloud_materials_DEFAULT.txt'
            print(f'{path}\n')
        else:
            print('uh oh, invalid file_str')
            quit()

    # if crucial file, also check that the file exists
    if file_str == 'infile' or file_str == 'longitudes':

        path = f'{loc}/{name}.txt'

        if os.path.isfile(path):
            print(f"The file exists at: {path}\n")
            return path
        else:
            print(f"The file does not exist at: {path}")
            quit()
    
    return path

def check_for_outpath(file_locs_and_names_dict,file_str):
    if f'{file_str}_loc' in file_locs_and_names_dict:
        loc = file_locs_and_names_dict[f'{file_str}_loc']

        if not os.path.exists(loc):
            os.makedirs(loc)
            print(f"Directory created at: {loc}")
        else:
            print(f"Directory already exists at: {loc}")

        print(f'Will output all plots and npz files to location: {loc}\n')

    else:
        print(f'No {file_str}_loc for read in provided,')
        print('will output all plots and npz files local to where the script is run (i.e. ./).')
        print('If this is not your desired effect,')
        print('please check your file_names_and_locs.txt file\n')
        loc = './'

    return loc

def check_for_run_name(file_locs_and_names_dict,file_str):
    if f'{file_str}' in file_locs_and_names_dict:
        name = file_locs_and_names_dict[f'{file_str}']
    else:
        print(f'No {file_str} provided,')
        print('Will use the infile_name as run_name\n')
        name = file_locs_and_names_dict['infile_name']

    return name