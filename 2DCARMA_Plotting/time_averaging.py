import numpy as np
from save_arrays import output_dictionary_to_compressed_npzfile
  
########################################################
# routine to time average the 2D read in arrays
def do_time_averaging(temporal_global_cloud_dict,mmr_chemical_element_dict,svp_chemical_element_dict,longitude_count_array,pressure_array,longitudes_array,r_array,infile_name):
    print('Starting time averaging')

    print('Longitude count array:',longitude_count_array)

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
    time_averaged_mmr_chemical_element_dict = {key: np.sum(value, axis=dim_to_sum)/longitude_count_array for key,value in mmr_chemical_element_dict.items()}
    time_averaged_svp_chemical_element_dict = {key: np.sum(value, axis=dim_to_sum)/longitude_count_array for key,value in svp_chemical_element_dict.items()}

    print('Time averaging complete')

    ########################################################
    # Save all the timeaveraged arrays here too, 
    # then don't have to pass huge number of arrays

    dicts_to_save_dict = {'time_averaged_global_cloud_properties' : time_averaged_global_cloud_dict,
                            'time_averaged_mmr_chemical_element_properties' : time_averaged_mmr_chemical_element_dict,
                            'time_averaged_svp_chemical_element_properties' : time_averaged_svp_chemical_element_dict,
                        }

    saved_time_averaged_dict_names_list = []
    for name, dict in dicts_to_save_dict.items():
        saved_dict_name = output_dictionary_to_compressed_npzfile(pressure_array,longitudes_array,r_array,longitude_count_array,
                                dict,name,infile_name)
        saved_time_averaged_dict_names_list.append(saved_dict_name)

    return saved_time_averaged_dict_names_list

    '''
    # Old way of doing it, saving until testing is done.
    time_averaged_global_cloud_dict = {}
    time_averaged_mmr_chemical_element_dict = {}
    time_averaged_svp_chemical_element_dict = {}
    for key, value in temporal_global_cloud_dict.items():
        print(np.shape(value))
        time_averaged_global_cloud_dict[key] = np.sum(value, axis=dim_to_sum)
        time_averaged_global_cloud_dict[key] /= longitude_count_array
        print(np.shape(time_averaged_global_cloud_dict[key]))
    for key, value in mmr_chemical_element_dict.items(): # These dictionaries share the same keys
        time_averaged_mmr_chemical_element_dict[key] = np.sum(value, axis=dim_to_sum)
        time_averaged_mmr_chemical_element_dict[key] /= longitude_count_array
        
        time_averaged_svp_chemical_element_dict[key] = np.sum(svp_chemical_element_dict[key], axis=dim_to_sum)
        time_averaged_svp_chemical_element_dict[key] /= longitude_count_array
    '''