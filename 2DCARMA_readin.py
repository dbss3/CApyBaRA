########################################################
# Import neccessary libraries, trying to keep this minimal as always
import numpy as np
import time

########################################################
# constants and setup values for plots
# assume these are material densities, not sure of the units
pressure_to_bar = 1e-6

M_TIO2 = 4.23
M_FE = 7.874
M_KCL = 1.98
M_ZNS = 4.1
M_NA2S = 1.856
M_MNS = 3.3
M_CR = 7.19
M_MG2SIO4 = 3.21
M_AL2O3 = 3.95

'''
plt.rcParams['figure.figsize'] = (12,6)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 20
'''

levels_num = np.array([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8])/1e7
levels_mass = np.array([1e-16,1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])/1e7

# For now we need to hard-code the names of these groups
# because there is no header, but in the future I will 
# add headers for this to 2DCARMA
group_name_list = ['tio2',
			'het_fe',
			'het_mg2sio4',
			'pure_fe',
			'het_cr',
			'cr',
			'mns',
			'na2s',
			'zns',
			'kcl',
			'al2o3'
			]

# This looks like Diana is just approximating that the cores << mantle
# for heterogeneous particles, as assuming they just have the mantle density.
# This is certainly true for larger grain sizes, but technically varies with bin size
# BUT, given the material densities are < factor 8 from each other, shouldn't matter much

group_density_list = [
M_TIO2, # 'tio2',
M_FE, # 'het_fe',
M_MG2SIO4, # 'het_mg2sio4',
M_FE, # 'pure_fe',
M_CR, # 'het_cr',
M_CR, # 'cr',
M_MNS, # 'mns',
M_NA2S, # 'na2s',
M_ZNS, # 'zns',
M_KCL, # 'kcl',
M_AL2O3, # 'al2o3'
]

# Different to Diana, I'm going to plot all the groups seperately,
# for now with the same colours for pure and heterogeneous
# This was Diana's list:
# t_tio2:'Blues'
# t_al2o3:'Greens'
# t_het_mg2sio4:'Purples'
# t_pure_fe+t_het_fe:'Reds' !!!!!
# t_cr:'Oranges'
# t_mns:'copper'
# t_na2s:'bone'
# t_kcl:'pink'
# t_zns:'Greys'

group_colours_list = [
'Blues',   # 'tio2',
'Reds',    # 'het_fe',
'Purples', # 'het_mg2sio4',
'Reds',    # 'pure_fe',
'Oranges', # 'het_cr',
'Oranges', # 'cr',
'copper',  # 'mns',
'bone',    # 'na2s',
'Greys',   # 'zns',
'pink',    # 'kcl',
'Greens'   # 'al2o3'
]

# Create group_properties_dict
groups_properties_dict = {name: {'group_density': density, 'group_colour': colour} for name, density, colour in zip(group_name_list, group_density_list, group_colours_list)}


#Also hoping to not have this hardcoded in the future
chemical_element_list = ['ti',
							'fe',
							'mg',
							'cr',
							'zn',
							'al',
							'kc'
						]

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

	# Now lets setup a dictionary for time-averaged 3D arrays of the above:
	time_averaged_global_cloud_array_3D = np.zeros((nz,nbin,ilong))

	time_averaged_global_cloud_dict = {}
	for g, group_name in enumerate(group_name_list):
		time_averaged_global_cloud_dict[group_name] = np.copy(time_averaged_global_cloud_array_3D)

	return temporal_global_cloud_dict, time_averaged_global_cloud_dict

def setup_dict_3D_and_time_averaged_dict_2D(nz,ilong,ntime):
	# Setup 3D arrays for pressure, longitude, time
	temporal_global_array_3D = np.zeros(((nz,ilong,ntime)))

	#believe these are for the mass mixing ratios of the gases, and the species vapour pressure
	mmr_chemical_element_dict = {}
	svp_chemical_element_dict = {}
	for c, chemical_element in enumerate(chemical_element_list): #Also see hard coding at top
		mmr_chemical_element_dict[chemical_element] = np.copy(temporal_global_array_3D)
		svp_chemical_element_dict[chemical_element] = np.copy(temporal_global_array_3D)

	# Now lets setup a dictionary for time-averaged 2D arrays of the above:
	time_averaged_chemical_element_array_2D = np.zeros((nz,ilong))

	time_averaged_mmr_chemical_element_dict = {}
	time_averaged_svp_chemical_element_dict = {}
	for c, chemical_element in enumerate(chemical_element_list):
		time_averaged_mmr_chemical_element_dict[chemical_element] = np.copy(time_averaged_chemical_element_array_2D)
		time_averaged_svp_chemical_element_dict[chemical_element] = np.copy(time_averaged_chemical_element_array_2D)

	return mmr_chemical_element_dict, svp_chemical_element_dict, time_averaged_mmr_chemical_element_dict, time_averaged_svp_chemical_element_dict

############################################################
# I have discovered that these arrays are not used, at least in this file
'''
def setup_dict_2D(nz, nbin):
	# For presssure and particle size
	pressure_and_particle_size_array_2D = np.zeros((nz,nbin))

	#Now for just chemical elements
	cloud_chemical_element_list = ['ti',
								   'fe',
								   'mg',
								   'cr',
								   'mn',
								   'na',
								   'zn',
								   'k',
								   'al',
								   'pure_cr',
								   'pure_fe'
								  ]

	# I'm not totally clear on what these are either,
	# hence the not so clear array names here
	p_chemical_element_dict = {}
	pm_chemical_element_dict = {}
	for c, cloud_chemical_element in enumerate(cloud_chemical_element_list):
		p_chemical_element_dict[cloud_chemical_element] = np.copy([pressure_and_particle_size_array_2D])
		pm_chemical_element_dict[cloud_chemical_element] = np.copy([pressure_and_particle_size_array_2D])

	return p_chemical_element_dict, pm_chemical_element_dict
'''

########################################################
# Functions to do the actual read in of the files
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
	temporal_global_cloud_dict, time_averaged_global_cloud_dict = setup_dict_4D_and_time_averaged_dict_3D(nz,nbin,ilong,ntime)
	mmr_chemical_element_dict, svp_chemical_element_dict, time_averaged_mmr_chemical_element_dict, time_averaged_svp_chemical_element_dict = setup_dict_3D_and_time_averaged_dict_2D(nz,ilong,ntime)

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
					temporal_global_cloud_dict[group_name][k,j,long_index,i] /= (math.log(ru_array[j,0]) - math.log(rl_array[j,0]))
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

	########################################################
	#TODO: space here to save the full data to compressed npz and then its like 0.1s to read in hopefully

	########################################################
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
		saved_dict_name = output_dictionary_to_compressed_npzfile(pressure_array,longitudes_array,r_array,
								dict,name)
		saved_time_averaged_dict_names_list.append(saved_dict_name)

	return saved_time_averaged_dict_names_list

########################################################
#function to output the data
def output_dictionary_to_compressed_npzfile(pressure_array,longitudes_array,r_array,data_dict,outfile_name):
	# Construct output dictionary
	output_dict = {
					'pressure_array':pressure_array,
					'longitudes_array':longitudes_array,
					'r_array':r_array,
					**data_dict # Unpack the dictionary into the new dictionary
	}

	saved_dict_name = f'{outfile_name}_{infile_name}.npz'

	np.savez_compressed(saved_dict_name, 
						**output_dict
						)
	return saved_dict_name