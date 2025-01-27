import numpy as np
from save_arrays import output_dictionary_to_compressed_npzfile
  
########################################################
# routine to time average the 2D read in arrays
def do_time_averaging(loaded_dicts_list,outfile_loc,run_name):
    print('Unpack the saved dict')
    # Need the axis arrays, i.e. pressure, longitude, long_count, r
    # for this just going to use the first one (should be cloud_properties)
    axis_arrays_dict = loaded_dicts_list[0]
    
    # Ok now unpack the axis arrays from this dictionary
    pressure_array = axis_arrays_dict['pressure_array']
    longitudes_array = axis_arrays_dict['longitudes_array']
    r_array = axis_arrays_dict['r_array']
    longitude_count_array = axis_arrays_dict['longitude_count_array']
    time_array = axis_arrays_dict['time_array']

    print('Longitude count array:',longitude_count_array)

    # Ok and then need to reconstruct the other dictionaries from the saved file
    temporal_global_cloud_dict = {key:value for key,value in loaded_dicts_list[0].items()  if '_array' not in key}
    mmr_cloud_element_dict     = {key:value for key,value in loaded_dicts_list[1].items()  if '_array' not in key}
    svp_cloud_element_dict     = {key:value for key,value in loaded_dicts_list[2].items()  if '_array' not in key}

    print('Finished unpacking saved npz files')

    print('Starting time averaging')
    # Now doing the time averaging
    # Dimension to sum over
    dim_to_sum = -1 # time is the last dimension, for both 3D and 4D temporal arrays, so sum over that

    # Create the new dictionaries with summed arrays
    # Doing explicit for loop, less repetition and clearer
    # As each longitude has a different number of saved steps, need to normalise
    # temporal_global_cloud_dict(nz,nbin,ilong,ntime), summing over last axis
    # time_averaged_global_cloud_dict(nz,nbin,ilong)

    #I hope this concise format works, my only concern is if the division is applied to the correct axis
    time_averaged_global_cloud_dict = {key: np.sum(value, axis=dim_to_sum)/longitude_count_array for key,value in temporal_global_cloud_dict.items()}
    time_averaged_mmr_cloud_element_dict = {key: np.sum(value, axis=dim_to_sum)/longitude_count_array for key,value in mmr_cloud_element_dict.items()}
    time_averaged_svp_cloud_element_dict = {key: np.sum(value, axis=dim_to_sum)/longitude_count_array for key,value in svp_cloud_element_dict.items()}

    print('Time averaging complete')

    ########################################################
    # Save all the timeaveraged arrays here too, 
    # then don't have to pass huge number of arrays

    print('Saving time averaged dicts')

    dicts_to_save_dict = {'time_averaged_global_cloud_properties' : time_averaged_global_cloud_dict,
                            'time_averaged_mmr_cloud_element_properties' : time_averaged_mmr_cloud_element_dict,
                            'time_averaged_svp_cloud_element_properties' : time_averaged_svp_cloud_element_dict,
                        }

    saved_time_averaged_dict_paths_list = []
    for name, dict in dicts_to_save_dict.items():
        saved_dict_path = output_dictionary_to_compressed_npzfile(pressure_array,longitudes_array,r_array,longitude_count_array,time_array,
                                dict,name,outfile_loc,run_name)
        saved_time_averaged_dict_paths_list.append(saved_dict_path)

    print('finished saving')

    return saved_time_averaged_dict_paths_list

    '''
    # Old way of doing it, saving until testing is done.
    time_averaged_global_cloud_dict = {}
    time_averaged_mmr_cloud_element_dict = {}
    time_averaged_svp_cloud_element_dict = {}
    for key, value in temporal_global_cloud_dict.items():
        print(np.shape(value))
        time_averaged_global_cloud_dict[key] = np.sum(value, axis=dim_to_sum)
        time_averaged_global_cloud_dict[key] /= longitude_count_array
        print(np.shape(time_averaged_global_cloud_dict[key]))
    for key, value in mmr_cloud_element_dict.items(): # These dictionaries share the same keys
        time_averaged_mmr_cloud_element_dict[key] = np.sum(value, axis=dim_to_sum)
        time_averaged_mmr_cloud_element_dict[key] /= longitude_count_array
        
        time_averaged_svp_cloud_element_dict[key] = np.sum(svp_cloud_element_dict[key], axis=dim_to_sum)
        time_averaged_svp_cloud_element_dict[key] /= longitude_count_array
    '''