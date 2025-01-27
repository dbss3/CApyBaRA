import numpy as np

########################################################
#function to output the arrays to a compressed npz file
def output_dictionary_to_compressed_npzfile(pressure_array,longitudes_array,r_array,longitude_count_array,data_dict,dict_name,outfile_loc,run_name):
	# Construct output dictionary
	output_dict = { 'pressure_array':pressure_array,
					'longitudes_array':longitudes_array,
					'r_array':r_array,
					'longitude_count_array':longitude_count_array,
					**data_dict # Unpack the dictionary into the new dictionary
	}

	saved_dict_path = f'{outfile_loc}/{dict_name}_{run_name}.npz'

	np.savez_compressed(saved_dict_path, 
						**output_dict
						)
	return saved_dict_path