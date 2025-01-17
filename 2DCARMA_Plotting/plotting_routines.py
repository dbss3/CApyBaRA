import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from constants import pressure_to_bar

########################################################
# Scripts that do the calculations for the various kinds of plots

def unpacking_arrays_from_loaded_dict(loaded_dict,group_properties_dict,group_name):
	#unpack some of the arrays from the dictionary
	pressure_array = loaded_dict['pressure_array']

	#presures are in dyn/cm^2 convert to bar
	pressure_array *= pressure_to_bar

	longitudes_array = loaded_dict['longitudes_array']
	r_array = loaded_dict['r_array']
	#print(pressure_array)

	#Need to determine axis for the group to select the appropriate bin sizes and the mass number density
	group_name_list = group_properties_dict.keys()
	group_index = group_name_list.index(group_name)

	time_averaged_global_specific_group_array = loaded_dict[group_name]
	specific_group_r_array = r_array[:,group_index]
	specific_group_properties_dict = group_properties_dict[group_name]
	print(specific_group_properties_dict.items())

	return time_averaged_global_specific_group_array, specific_group_properties_dict, specific_group_r_array, pressure_array, longitudes_array

def handle_zeros_and_make_levels(contour_value_array):

	# Also handle zeros - contour can't handle 0 in log
	contour_value_array[contour_value_array==0] = np.nan

	is_all_nan = np.all(np.isnan(contour_value_array))
	if np.all(is_all_nan):
		print('None of this cloud species, forcing blank plot')
		# Replace the nans back with 1e-99 because contour can't handle all nan
		# Dummy levels
		levels = np.array([1,2])
	else:
		print('There is some of this cloud species! Plotting :)')
		# Determine the min and max
		max_value = np.nanmax(contour_value_array)
		min_value = np.nanmin(contour_value_array)
		#print(contour_value_array)
		print(f'max and min of contour: {max_value} {min_value}')
		if min_value < 1e-20:
			min_value = 1e-20
		# make levels array:
		ceil_max_log_value = int(np.ceil(np.log10(max_value)))
		flor_min_log_value = int(np.floor(np.log10(min_value)))
		int_diff = int(ceil_max_log_value-flor_min_log_value)
		print(flor_min_log_value,ceil_max_log_value,int_diff)

		levels = np.logspace(flor_min_log_value,ceil_max_log_value,int_diff+1)
		print(levels)
	
	return contour_value_array, levels

def cloud_contour_plot_cartesian_longitude_pressure(axis,group_name,loaded_dict,group_properties_dict,mass_or_numb_density = 'numb',alpha=0.75):
	# plan is to make a function that handles the plotting of the contours,
	# and hand it the axis to plot to, so can stack materials easily and also
	# make their own plots by different calls to this function

	#unpack some of the arrays from the dictionary
	time_averaged_global_specific_group_array, specific_group_properties_dict,\
	specific_group_r_array, pressure_array, longitudes_array = \
	unpacking_arrays_from_loaded_dict(loaded_dict,group_properties_dict,group_name)

	# Gotta sum over the number density bins
	# So time_averaged_global_specific_group_array has dimensions (nz,nbin,ilong)
	# Need to sum up over the nbin dimension
	time_averaged_global_specific_group_summed_nbin_array = np.sum(time_averaged_global_specific_group_array, axis=1)
	# So this final array should have dimensions (nz,ilong)

	#renaming for convenience, when doing mass and number density need same name
	contour_value_array = time_averaged_global_specific_group_summed_nbin_array

	# handle zero values, which cause problems with log scale,
	# then also generate appropriate levels
	contour_value_array, levels = handle_zeros_and_make_levels(contour_value_array)

	# Construct X-Y grids for contourf
	X_grid, Y_grid = np.meshgrid(longitudes_array,pressure_array) 

	# Plot to the desired axis
	contour_object = axis.contourf(X_grid, Y_grid, contour_value_array,
				norm=LogNorm(), cmap=specific_group_properties_dict['group_colour'], 
				alpha=alpha, levels = levels)

	return contour_object, X_grid, Y_grid, contour_value_array

def cloud_contour_plot_cartesian_longitude_particle_bin(axis,group_name,desired_pressure,loaded_dict,group_properties_dict,mass_or_numb_density = 'numb',alpha=0.75):
	# plan is to make a function that handles the plotting of the contours,
	# and hand it the axis to plot to, so can stack materials easily and also
	# make their own plots by different calls to this function

	#unpack some of the arrays from the dictionary
	time_averaged_global_specific_group_array, specific_group_properties_dict,\
	specific_group_r_array, pressure_array, longitudes_array = \
	unpacking_arrays_from_loaded_dict(loaded_dict,group_properties_dict,group_name)

	#And now need to pick the level closest to the desired pressure
	press_index = np.argmin(np.abs(pressure_array-desired_pressure)) #make sure they have the same units!

	time_averaged_specific_level_specific_group_array = time_averaged_global_specific_group_array[press_index]
	#An explanation: This last array has 2 dimensions: (nbin,ilong)

	# Now determine the particle_mass_density array OR *_number_denisty_array
	# For the selected group index, and pressure_index
	# choice of mass or number density
	if mass_or_numb_density == 'numb':
		contour_value_array = time_averaged_specific_level_specific_group_array
	elif mass_or_numb_density == 'mass':
		# Have to convert using material_density*4pi/3*r^3
		contour_value_array = time_averaged_specific_level_specific_group_array*(specific_group_properties_dict['group_density']*((4./3.)*np.pi*(specific_group_r_array*1e-4)**3.))
	else:
		print('Invalid value for mass_or_numb_density, please pick mass or numb, not: ',mass_or_numb_density)
	
	# handle zero values, which cause problems with log scale,
	# then also generate appropriate levels
	contour_value_array, levels = handle_zeros_and_make_levels(contour_value_array)

	# Construct X-Y grids for contourf
	X_grid, Y_grid = np.meshgrid(longitudes_array,specific_group_r_array) #the index here is the bin index and then group index
	
	# Plot to the desired axis
	contour_object = axis.contourf(X_grid, Y_grid, contour_value_array,
			      norm=LogNorm(), cmap=specific_group_properties_dict['group_colour'], 
				  alpha=alpha) #, levels = levels)

	return contour_object, X_grid, Y_grid, contour_value_array

def cloud_contour_plot_cartesian_particle_bin_pressure(axis,group_name,desired_longitude,loaded_dict,group_properties_dict,mass_or_numb_density = 'numb',alpha=0.75,opening_angle=20):
	# plan is to make a function that handles the plotting of the contours,
	# and hand it the axis to plot to, so can stack materials easily and also
	# make their own plots by different calls to this function

	#unpack some of the arrays from the dictionary
	time_averaged_global_specific_group_array, specific_group_properties_dict,\
	specific_group_r_array, pressure_array, longitudes_array = \
	unpacking_arrays_from_loaded_dict(loaded_dict,group_properties_dict,group_name)

	# time_averaged_global_specific_group_array dimensions: (nz,nbin,ilong)
	# And now masks +- opening_angle around the desired_longitude
	# Indexing on the longitude_array values
	mask = np.where((longitudes_array < desired_longitude+opening_angle) & (longitudes_array > desired_longitude-opening_angle))[0]
	print(longitudes_array)
	print(mask)

	longitude_specific_masked_array = time_averaged_global_specific_group_array[:,:,mask]

	# OK, now we want to sum over the remaining values of the longitude dimension
	contour_value_array = np.sum(longitude_specific_masked_array, axis=-1)

	# handle zero values, which cause problems with log scale,
	# then also generate appropriate levels
	contour_value_array, levels = handle_zeros_and_make_levels(contour_value_array)

	# Construct X-Y grids for contourf
	X_grid, Y_grid = np.meshgrid(specific_group_r_array,pressure_array) #the index here is the bin index and then group index
	
	# Plot to the desired axis
	contour_object = axis.contourf(X_grid, Y_grid, contour_value_array,
			      norm=LogNorm(), cmap=specific_group_properties_dict['group_colour'], 
				  alpha=alpha) #, levels = levels)

	return contour_object, X_grid, Y_grid, contour_value_array

########################################################
# plotting wrapping routines
def plotter(group_name,desire_value,loaded_dict,which_plot):
	# save the figures, WANT TO also make it so that this
	# saves the arrays that are plotted, so axes etc can be 
	# fixed more easily
	fig, ax = plt.subplots(nrows = 1, ncols = 1)

	if which_plot == 'longitude_particle_bin_contour':
		print(f'Plotting lon vs num den for {group_name}')
		contour_object, X_grid, Y_grid, contour_value_array = \
			cloud_contour_plot_cartesian_longitude_particle_bin(ax,group_name,desire_value,loaded_dict)
		cbar = plt.colorbar(contour_object)
		cbar.set_label(r'Number Density [cm$^{-3}$]', rotation=-90, labelpad=30) #fontsize=20
		ax.set_xlabel(r'Longitude [degree]')
		ax.set_ylabel(r'Particle Radius [$\mu$m]')
		ax.set_yscale('log')
		#ax.set_xlim(min(longitudes_array),max(longitudes_array))
		#ax.set_ylim(1,min(pressure_array/1.e6))
		save_name = f'particle_bin_of_{group_name}_at_{desire_value}bar'


	elif which_plot == 'longitude_pressure_contour':
		print(f'Plotting lon vs pressure for {group_name}')
		contour_object, X_grid, Y_grid, contour_value_array = \
			cloud_contour_plot_cartesian_longitude_pressure(ax,group_name,loaded_dict)
		cbar = plt.colorbar(contour_object)
		cbar.set_label(r'Integrated Number Density [--]', rotation=-90, labelpad=30) #fontsize=20
		ax.set_xlabel(r'Longitude [degree]')
		ax.set_ylabel(r'Pressure [bar]')
		ax.set_yscale('log')
		#ax.set_xlim(min(longitudes_array),max(longitudes_array))
		#ax.set_ylim(1,min(Y_grid/1.e6))
		ax.axvline(90,ls='--',color='g')
		ax.axvline(-90,ls='--',color='g')
		lon_ticks = np.arange(180,-180,-45) # doing backwards because lons are 180 >= l > -180
		ax.set_xticks(lon_ticks)
		ax.set_xticklabels(lon_ticks)
		plt.gca().invert_yaxis()
		save_name = f'integrated_number_density_of_{group_name}'

	elif which_plot == 'particle_bin_pressure_contour':
		print(f'Plotting num den vs press for {group_name}')
		contour_object, X_grid, Y_grid, contour_value_array = \
			cloud_contour_plot_cartesian_particle_bin_pressure(ax,group_name,desire_value,loaded_dict)
		cbar = plt.colorbar(contour_object)
		cbar.set_label(r'Integrated Number Density [--]', rotation=-90, labelpad=30) #fontsize=20
		ax.set_xlabel(r'Particle Radius [$\mu$m]')
		ax.set_ylabel(r'Pressure [bar]')
		ax.set_yscale('log')
		ax.set_xscale('log')
		#ax.set_xlim(min(longitudes_array),max(longitudes_array))
		#ax.set_ylim(1,min(Y_grid/1.e6))
		plt.gca().invert_yaxis()
		save_name = f'integrated_number_density_of_{group_name}_at_{desire_value}deg'
	else:
		print('Please pick a valid type of plot:')
		print('')
	
	# Save the figures
	plt.savefig(f'{save_name}.pdf', format = 'pdf', bbox_inches='tight')
	print('Figure Saved \n')

	# Close whatever the figure plotted was
	plt.close()

	# Save the arrays:
	np.savez_compressed(save_name,
						X_grid = X_grid,
						Y_grid = Y_grid,
						contour_value_array = contour_value_array
						)

	return