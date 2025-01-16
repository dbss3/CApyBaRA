import numpy as np

########################################################
# Functions to do setup of the arrays

def setup_arrays_1D(nz,ntime,ilong):
	# Setup 1D arrays 

	# for the p-T-z structure as 1D arrays
	profile_array_1D = np.zeros(nz)

	pressure_array = np.copy(profile_array_1D) # p in Diana's code
	temperature_array = np.copy(profile_array_1D) # t in Diana's code
	height_array = np.copy(profile_array_1D) # z in Diana's code

	# for the temporal variation with rotation
	temporal_array_1D = np.zeros(ntime)

	time_array     = np.copy(temporal_array_1D)
	distance_array = np.copy(temporal_array_1D)
	rotation_array = np.copy(temporal_array_1D)

	# for longitude summation
	longitude_count_array = np.zeros((ilong))

	return pressure_array, temperature_array, height_array, time_array, distance_array, rotation_array, longitude_count_array

def setup_arrays_2D(nbin,ngroup):
	# Setup 2D arrays 

	# for the different groups
	particle_size_and_composition_array_2D = np.zeros((nbin,ngroup))

	# Not entirely clear what these are yet, I'm sure it'll become
	# clear as I get further in to Diana's script.
	# I'm guessing r is radius, rl is radius_lower, ru is radius_upper
	# ms is mass, dr is the bin width in r?
	r_array  = np.copy(particle_size_and_composition_array_2D)
	ms_array = np.copy(particle_size_and_composition_array_2D)
	dr_array = np.copy(particle_size_and_composition_array_2D)
	rl_array = np.copy(particle_size_and_composition_array_2D)
	ru_array = np.copy(particle_size_and_composition_array_2D)

	return r_array , ms_array, dr_array, rl_array, ru_array

def setup_dict_4D_and_time_averaged_dict_3D(nz,nbin,ilong,ntime):
	# Setup 4D arrays for pressure, particle size, long, and time
	temporal_global_cloud_array_4D = np.zeros((nz,nbin,ilong,ntime))
	
	temporal_global_cloud_dict = {}
	for g, group_name in enumerate(group_name_list): #See hardcoding at the top
		temporal_global_cloud_dict[group_name] = np.copy(temporal_global_cloud_array_4D)
		
    '''
	# Getting rid of time-averaging here, as I don't think its needed
	# But keeping the code until I've tested this.
	# Now lets setup a dictionary for time-averaged 3D arrays of the above:
	time_averaged_global_cloud_array_3D = np.zeros((nz,nbin,ilong))

	time_averaged_global_cloud_dict = {}
	for g, group_name in enumerate(group_name_list):
		time_averaged_global_cloud_dict[group_name] = np.copy(time_averaged_global_cloud_array_3D)
    '''

	return temporal_global_cloud_dict #, time_averaged_global_cloud_dict

def setup_dict_3D_and_time_averaged_dict_2D(nz,ilong,ntime):
	# Setup 3D arrays for pressure, longitude, time
	temporal_global_array_3D = np.zeros(((nz,ilong,ntime)))

	#believe these are for the mass mixing ratios of the gases, and the species vapour pressure
	mmr_chemical_element_dict = {}
	svp_chemical_element_dict = {}
	for c, chemical_element in enumerate(chemical_element_list): #Also see hard coding at top
		mmr_chemical_element_dict[chemical_element] = np.copy(temporal_global_array_3D)
		svp_chemical_element_dict[chemical_element] = np.copy(temporal_global_array_3D)

	'''
	# Getting rid of time-averaging here, as I don't think its needed
	# But keeping the code until I've tested this.
	# Now lets setup a dictionary for time-averaged 2D arrays of the above:
	time_averaged_chemical_element_array_2D = np.zeros((nz,ilong))

	time_averaged_mmr_chemical_element_dict = {}
	time_averaged_svp_chemical_element_dict = {}
	for c, chemical_element in enumerate(chemical_element_list):
		time_averaged_mmr_chemical_element_dict[chemical_element] = np.copy(time_averaged_chemical_element_array_2D)
		time_averaged_svp_chemical_element_dict[chemical_element] = np.copy(time_averaged_chemical_element_array_2D)
		'''

	return mmr_chemical_element_dict, svp_chemical_element_dict #, time_averaged_mmr_chemical_element_dict, time_averaged_svp_chemical_element_dict