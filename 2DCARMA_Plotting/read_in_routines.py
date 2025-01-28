import numpy as np
from setup_arrays_routines import setup_arrays_1D, setup_arrays_2D,setup_dicts_3D,setup_dicts_4D
from constants import pressure_to_bar
from save_arrays import output_dictionary_to_compressed_npzfile

########################################################
# Functions to do the actual read in of the files

def read_file_to_dict(file_path):
	"""
	Reads a file and parses it into a dictionary. Lines starting with '#' are ignored,
	and lines starting with '!' are treated as keys, with the following line as the corresponding value.

	:param file_path: Path to the input file
	:return: Dictionary with keys and values parsed from the file
	"""
	file_locs_and_names_dict = {}
	
	with open(file_path, 'r') as file:
		lines = file.readlines()
		
	# Clean up lines and filter out comments
	clean_lines = []
	for line in lines:
		stripped_line = line.split('#')[0].strip()  # Remove comments and trim whitespace
		if stripped_line:  # Skip empty lines
			clean_lines.append(stripped_line)
	
	# Iterate through lines to build the dictionary
	for clean_line in clean_lines:
		if clean_line.startswith('!'):
			key = clean_line[1:].strip()  # Remove '!' and strip whitespace
			file_locs_and_names_dict[key] = []
		else: #means its a value line, hopefully the previous cleaning has worked
			value = clean_line.strip()  # Strip whitespace
			file_locs_and_names_dict[key].append(value)
	
	#Now go through the dict and clean up the lists to be strings unless multi values
	for key, value in file_locs_and_names_dict.items():
		if len(value) == 0:
			print(f'Error: key {key} has no specified value')
			quit()
		elif len(value) == 1: #only one entry so setting to string
			print(value[0])
			file_locs_and_names_dict[key] = value[0]
		else:
			#leaving list
			pass
	
	return file_locs_and_names_dict

def load_txt_file(file_name):
    """
    Reads a text file and returns a list of strings from each line,
    ignoring lines that start with '#'.
	# Generated using CHAT GPT
    
    Parameters:
        file_name (str): The name of the text file to read.
    
    Returns:
        list: A list of strings for each line not starting with '#'.
    """
    result = []
    try:
        with open(file_name, 'r') as file:
            for line in file:
                stripped_line = line.strip()
                if not stripped_line.startswith("#"):
                    result.append(stripped_line)
    except FileNotFoundError:
        print(f"Error: The file {file_name} was not found.")
    return result


def read_cloud_properties(file_name):
	"""
	Parses a structured text file into a dictionary based on a key header line,
	ignoring comments (#) and processing the cloud_elements_list as a list of strings.

	Parameters:
		file_name (str): The name of the text file to parse.

	Returns:
		dict: A dictionary where keys are taken from the header line starting with '!',
				and values are lists of strings parsed from the lines.
	"""
	header_keys = None
	data_list = []

	try:
		with open(file_name, 'r') as file:
			for line in file:
				# Remove comments starting with #
				line = line.split("#")[0].strip()
				
				# Ignore empty lines
				if not line:
					continue
				
				# Identify header key
				if line.startswith("!"):
					header_keys = line[1:].strip().split(',')  # Remove '!' and strip whitespace
					continue #continue means loop not continue, stupid
				
				# Process data lines
				if len(header_keys) != 0:
					parts = line.split(",")
					# Strip whitespace from each part
					parts = [part.strip() for part in parts]
					
					# Extract and process the cloud_elements_list (last part)
					if parts[-1].startswith("[") and parts[-1].endswith("]"):
						cloud_elements = parts[-1][1:-1].split(",")  # Remove brackets and split
						cloud_elements = [el.strip() for el in cloud_elements]  # Strip each element
						parts[-1] = cloud_elements
					
					data_list.append(parts)
					#print(parts)

			# finally combine into an actual dictionary
			print(header_keys)
			print(data_list)
			data_list_transpose = [list(row) for row in zip(*data_list)]
			print(data_list_transpose)
			#have to transpose the list somehow (can't use numpy)
			#to get this
			parsed_data_dict = dict(zip(header_keys,data_list_transpose))
			for key, value in parsed_data_dict.items():
				print(key)
				print(value)
				print('')

	except FileNotFoundError:
		print(f"Error: The file {file_name} was not found.")
	except Exception as e:
		print(f"An error occurred: {e}")

	return parsed_data_dict

def read_in_for_2DCARMA(infile_path,longitudes_path,outfile_loc,run_name,cloud_properties_file_path,cloud_materials_file_path):
	########################################################
	# For now I need the to read in the longitudes seperately to get 
	# the number of longitudes, I will be adding this to 2DCARMA output
	longitudes_array = np.loadtxt(longitudes_path)
	ilong = len(longitudes_array)

	########################################################
	# Have now set this up so that I need to read these files in here
	# to setup group names and materials in the cloud elements
	cloud_properties_dict = read_cloud_properties(cloud_properties_file_path)

	# The groups are their own individual entries in the dictionary
	group_name_list = cloud_properties_dict['group_name']
	print(f'Read in expected group names: {group_name_list}')
    
	# Have to do a bit more work to get the unique materials comprising the clouds
	# Will work in the future when the output has headers

	#But ordering matters here, for now just reading in an ordered seperate file,
	cloud_element_list = load_txt_file(cloud_materials_file_path)
    #TODO: OK THIS FILE IS REALLY CONFUSING, its just read in now for consistency
    # with the code, but these arrays I don't use atm, and I don't know what they are

	########################################################
	# The file this one reads in is the main one, 
	# not the *_temp*.txt files, or the *_rates.txt, or *_flux.txt files
	infile = open(infile_path,'r')

	########################################################
	# Just read the first line for the run details
	# nz,ngroup,nelem,nbin,ngas,nstep,iskip
	# nz = number of pressure/height levels
	# ngroup = number of groups (i.e. the full core+shell particles)
	# nelem = number of elements NOT chemical elements, 
	# 		  but the number of material components of the particles
	#		  i.e. Fe core, Fe mantle, TiO2 core, MgSiO3 mantle
	# nbin = number of size (mass) bins for the cloud particles
	# ngas = number of gas species included
	# nstep = TOTAL number of steps the run was set up for
	#		  NB: THERE IS A PROBLEM USING THIS IF THE RUN DIDN'T COMPLETE
	# iskip = number of steps which are skipped between outputs
	
	line = infile.readline().split()
	nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
      
	# Do a quick check to see if ngroup matches the cloud dict loaded
	print('Checking input file matches expected groups')
	print(f'ngroup from input file: {ngroup}, group properties for: {len(group_name_list)}')
	if len(group_name_list) == ngroup:
		print('File matches expected cloud group setup\n')
	else:
		print('File DOES NOT match expected cloud groups, please input a correct cloud group file,')
		print(f'cloud_group file used atm: {cloud_properties_file_path}\n')
		quit()

	# Do a little computation, NB: assumes the run completes
	# ntime = number of timesteps recorded
	ntime=int(nstep/iskip)
	print(nstep,iskip,ntime)

	print('Setting up the arrays')
	########################################################
	# Set up all the arrays
	pressure_array, temperature_array, height_array, \
		time_array, distance_array, rotation_array, longitude_count_array = setup_arrays_1D(nz,ntime,ilong)
	r_array , ms_array, dr_array, rl_array, ru_array = setup_arrays_2D(nbin,ngroup)
	#p_cloud_element_dict, pm_cloud_element_dict = setup_dict_2D(nz, nbin)
	temporal_global_cloud_dict = setup_dicts_4D(nz,nbin,ilong,ntime,group_name_list)
	temporal_mmr_cloud_element_dict, temporal_svp_cloud_element_dict = setup_dicts_3D(nz,ilong,ntime,cloud_element_list)

	print('Arrays set up')

	print('Reading the particle size bins')
	########################################################
	# First block in the file is the bin and group array details
	for i in range(ngroup):
		for j in range(nbin):
			line = infile.readline().split()
			r_array[j,i] = float(line[2])
			ms_array[j,i] = float(line[3])
			dr_array[j,i] = float(line[4])
			rl_array[j,i] = float(line[5])
			ru_array[j,i] = float(line[6])

	print('finished')

	print('Reading p-T-z structure')
	########################################################
	# Next block is the vertical p-T-z structure
	for i in range(nz):
		line = infile.readline().split()
		height_array[i] = float(line[1])
		pressure_array[i] = float(line[3])
		#_ = float(line[2]) # Not sure what this is, Diana doesn't use it: Check CARMA
		temperature_array[i] = float(line[4])
	
	print('finished')

	print('Now reading the individual timesteps')
	########################################################
	# Now we get onto blocks that are for each time-step
	# Each time there is a line which describes the time-step first
	for i in range(ntime):
		end_time_step = i
		timestep_details = infile.readline().split()
		print(timestep_details)

		# If there are no more timesteps, when run doesn't complete,
		# means this line will not run
		if len(timestep_details) == 0:
			print('looks like this run did not complete,')
			print(f'ending read-in at i = {end_time_step}, expected: {ntime}\n')
			break

		# If here, continuing with read in
		if i == 0:
			# For some reason the first one just has the first time-step time
			# If the file is from an OG run, this will just be 0
			# If its a restart it should be whatever that start time is?
			time_array[i] = float(timestep_details[0])
			long_index = 0
		else:
			# Otherwise there are three values, for the timestep.
			# The time, the distance?, and the rotation?
			# So not each timestep is at the same longitude!!!
			time_array[i] = float(timestep_details[0])
			distance_array[i] = float(timestep_details[1])
			rotation_array[i] = float(timestep_details[2])
			current_longitude = distance_array[i] - (len(longitudes_array)*rotation_array[i])
			if int(round(current_longitude,0)) == ilong: #changed this to ilong to avoid hardcoding - might need a minus 1 etc
				long_index = 0
			else:
				long_index = int(round(current_longitude,0))
		
		print('iteration: ',i)
		print('reading time:',time_array[i])
		print('longitude:',long_index)
		print('')

		# Saving number of instances of each longitude to an array for time average calc
		longitude_count_array[long_index] +=1
		
		# Finally for this timestep block we can actually read in the properites
		for j in range(nbin):
			for k in range(nz):
				#print(i,j,k)
				line = infile.readline().split()

				# I'm going to split the line into the two sections, the cloud and then the gas
				# Ignore the first two vlaues as these are just the bin number and the height
				starting = 2
				cloud_properties_for_line_array = line[starting:ngroup+starting]
				#print(cloud_properties_for_line_array)
				# The gas values are interleafed so stepping every other value here 
				mmr_for_line_array = line[ngroup+starting::2]
				svp_for_line_array = line[ngroup+starting+1::2]

				#lets check the line lengths align with our hardcoded expectations
				#print('Are lines the right length for the hardcoded values:')
				#print(len(group_name_list),len(cloud_properties_for_line_array))
				#print(len(cloud_element_list),len(mmr_for_line_array),len(svp_for_line_array))

				for g, group_name in enumerate(group_name_list):
					temporal_global_cloud_dict[group_name][k,j,long_index,i] = float(cloud_properties_for_line_array[g])
					#need to normalise or something
					temporal_global_cloud_dict[group_name][k,j,long_index,i] /= (np.log(ru_array[j,0]) - np.log(rl_array[j,0]))
					#Lets compute this later - NOW DONE JUST BEFORE PLOTTING
					#mass_tio2[k,j,long_index,i] = tio2[k,j,long_index,i]*(M_TIO2*((4./3.)*np.pi*(r[j,0]*1e-4)**3.))

				for c, cloud_element in enumerate(cloud_element_list):
					temporal_mmr_cloud_element_dict[cloud_element][k,long_index,i] = float(mmr_for_line_array[c])
					temporal_svp_cloud_element_dict[cloud_element][k,long_index,i] = float(svp_for_line_array[c])
					
					#need to covert: partial pressure
					#and convert to bar, have to do square because both values above are in the wrong units
					temporal_mmr_cloud_element_dict[cloud_element][k,long_index,i] *= pressure_array[k]*pressure_to_bar**2
					temporal_svp_cloud_element_dict[cloud_element][k,long_index,i] *= pressure_array[k]*pressure_to_bar**2

	print('Finally Read-in complete\n')

	print(time_array)

	# Need to trim down all the time arrays to the correct size,
	# Doesn't affect time averaging, but it will make combining them neater
	if end_time_step != ntime:
		print('Triming arrays as incomplete run, endtime = ', end_time_step)
		print('Total number of read in times: ',np.sum(longitude_count_array))
		print('Expected number of times to read in: ', ntime)

		# 1D Arrays that need to be trimmed
		print('Check array before: ',time_array)
		time_array     = time_array[:end_time_step]
		distance_array = distance_array[:end_time_step]
		rotation_array = rotation_array[:end_time_step]
		print('Check array after: ',time_array)

		# Dictionaries that need their arrays trimmed
		for key, array in temporal_global_cloud_dict.items():
			trimmed_array = array[:,:,:,:end_time_step]
			temporal_global_cloud_dict[key] = trimmed_array
		
		for key, array in temporal_mmr_cloud_element_dict.items():
			trimmed_array = array[:,:,:end_time_step]
			temporal_mmr_cloud_element_dict[key] = trimmed_array
			
		for key, array in temporal_svp_cloud_element_dict.items():
			trimmed_array = array[:,:,:end_time_step]
			temporal_svp_cloud_element_dict[key] = trimmed_array

		print('Triming complete')

      
	########################################################
	# Save the full data to compressed npz

	print('Begining saving data')

	dicts_to_save_dict = {'temporal_global_cloud_properties': temporal_global_cloud_dict,
					      'temporal_mmr_cloud_element_properties':  temporal_mmr_cloud_element_dict,
						  'temporal_svp_cloud_element_properties':  temporal_svp_cloud_element_dict,
					   }
	
	saved_dict_paths_list = []
	for dict_name, dict in dicts_to_save_dict.items():
		print(f'saving: {dict_name}')
		saved_dict_path = output_dictionary_to_compressed_npzfile(pressure_array,longitudes_array,r_array,longitude_count_array,time_array,
								dict,dict_name,outfile_loc,run_name)
		saved_dict_paths_list.append(saved_dict_path)
		print('saved\n')
	
	print('saving complete. Read in Complete\n')
	return saved_dict_paths_list