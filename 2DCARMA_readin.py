import numpy as np
from array_setup_routines import setup_arrays_1D, setup_arrays_2D,setup_dicts_3D,setup_dicts_4D

########################################################
# Functions to do the actual read in of the files

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


def read_in(file_loc,longitudes_loc):
	# For now I need the to read in the longitudes seperately to get 
	# the number of longitudes, I will be adding this to 2DCARMA output
	longitudes_array = np.loadtxt(longitudes_loc)
	ilong = len(longitudes_array)

	########################################################
	# The file this one reads in is the main one, 
	# not the *_temp*.txt files, or the *_rates.txt, or *_flux.txt files
	infile = open(file_loc,'r')

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

	# Do a little computation, NB: assumes the run completes
	# ntime = number of timesteps recorded
	ntime=int(nstep/iskip)
	print(nstep,iskip,ntime)

	print('Setting up the arrays')
	########################################################
	# Set up all the arrays
	pressure_array, temperature_array, height_array, time_array, distance_array, rotation_array, longitude_count_array = setup_arrays_1D(nz,ntime,ilong)
	r_array , ms_array, dr_array, rl_array, ru_array = setup_arrays_2D(nbin,ngroup)
	#p_chemical_element_dict, pm_chemical_element_dict = setup_dict_2D(nz, nbin)
	temporal_global_cloud_dict = setup_dicts_4D(nz,nbin,ilong,ntime,group_name_list_file_loc)
	mmr_chemical_element_dict, svp_chemical_element_dict = setup_dicts_3D(nz,ilong,ntime,chemical_element_list_file_loc)

	print('arrays set up')

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
	for i in range(ntime-1):
		if i == 0:
			# For some reason the first one just has the first time-step time
			# If the file is from an OG run, this will just be 0
			# If its a restart it should be whatever that start time is?
			time_array[i] = float(infile.readline())
			long_index = 0
		else:
			# Otherwise there are three values, for the timestep.
			# The time, the distance?, and the rotation?
			# So not each timestep is at the same longitude!!!
			timestep_details  = infile.readline().split()
			time_array[i] = float(timestep_details[0])
			distance_array[i] = float(timestep_details[1])
			rotation_array[i] = float(timestep_details[2])
			current_longitude = distance_array[i] - (len(longitudes_array)*rotation_array[i])
			if int(round(current_longitude,0)) == ilong: #changed this to ilong to avoid hardcoding - might need a minus 1 etc
				long_index = 0
			else:
				long_index = int(round(current_longitude,0))

		print('reading time:',time_array[i])
		print('longitude:',long_index)

		# Saving number of instances of each longitude to an array for time average calc
		longitude_count_array[long_index] +=1
		
		# Finally for this timestep block we can actually read in the properites
		for j in range(nbin):
			for k in range(nz):
				print(i,j,k)
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
				#print(len(chemical_element_list),len(mmr_for_line_array),len(svp_for_line_array))

				for g, group_name in enumerate(group_name_list):
					temporal_global_cloud_dict[group_name][k,j,long_index,i] = float(cloud_properties_for_line_array[g])
					#need to normalise or something
					temporal_global_cloud_dict[group_name][k,j,long_index,i] /= (np.log(ru_array[j,0]) - np.log(rl_array[j,0]))
					#Lets compute this later - NOW DONE JUST BEFORE PLOTTING
					#mass_tio2[k,j,long_index,i] = tio2[k,j,long_index,i]*(M_TIO2*((4./3.)*np.pi*(r[j,0]*1e-4)**3.))

				for c, chemical_element in enumerate(chemical_element_list):
					mmr_chemical_element_dict[chemical_element][k,long_index,i] = float(mmr_for_line_array[c])
					svp_chemical_element_dict[chemical_element][k,long_index,i] = float(svp_for_line_array[c])
					
					#need to covert: partial pressure
					#and convert to bar, have to do square because both values above are in the wrong units
					mmr_chemical_element_dict[chemical_element][k,long_index,i] *= pressure_array[k]*pressure_to_bar**2
					svp_chemical_element_dict[chemical_element][k,long_index,i] *= pressure_array[k]*pressure_to_bar**2

	print('Finally Read-in complete')