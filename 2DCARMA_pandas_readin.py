####################################################
import pandas as pd
import numpy as np
####################################################
file_path = './WASP-52b_clear_1xsolar_1Kzz.txt'

group_name_list = ['tio2',
            'al2o3',
			'het_fe',
			'het_mg2sio4',
			'het_cr',
			'mns',
			'na2s',
			'pure_fe',
			'cr',
			'kcl',
			'zns',
			]

#Also hoping to not have this hardcoded in the future
chemical_element_list = ['ti',
							'fe',
							'mg',
							'cr',
                            'dum1',
                            'dum2',
                            'dum3',
                            'kc',
							'al'
						]

M_dict = {'M_kcl' : 1.98,
          'M_zns' : 4.1,
            'M_na2s' : 1.856,
            'M_mns' : 3.3,
            'M_cr' : 7.19,
            'M_het_cr' : 7.19,
            'M_het_fe' : 7.874,
            'M_pure_fe' : 7.874,
            'M_het_mg2sio4' : 3.21,
            'M_tio2' : 4.23,
            'M_al2o3' : 3.95
}

gas_data_name_list = []
for c, chemical_element in enumerate(chemical_element_list):
    gas_data_name_list.append('mmr_%s'%chemical_element)
    gas_data_name_list.append('svp_%s'%chemical_element)

data_name_list = ['jbin','klevel']+group_name_list+gas_data_name_list

print(data_name_list)

lines_read = 0

longitude = np.loadtxt('../run/carma/longitudes_WASP-52b_clear_1xsolar.txt')
ilong = len(longitude)

####################################################
# Read just the first line of the file
with open(file_path, 'r') as file:
    first_line = file.readline().split()

# Parse the first line to extract the number of rows and columns
nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int, first_line)  # Adjust separator as per your data

print(nz,ngroup,nelem,nbin,ngas,nstep,iskip)

# Do a little computation, NB: assumes the run completes
# ntime = number of timesteps recorded
ntime=int(nstep/iskip)
print(nstep,iskip,ntime)

lines_read += 1

####################################################
# Define the number of lines to read
lines_to_read = ngroup * nbin  # Calculate the total number of lines to read

# Read the data into a Pandas DataFrame
particle_bin_df = pd.read_csv(
    file_path,
    sep='\s+',  # Automatically handles varying spaces as delimiters
    skiprows=lines_read,  # Skip the first line (metadata already read)
    nrows=lines_to_read,  # Read only the required number of lines
    header=None,  # No header in the file
    names=['igroup', 'jbin', 'r', 'ms', 'dr', 'rl', 'ru']  # Assign column names
)

for g, group_name in enumerate(group_name_list):
   particle_bin_df.loc[particle_bin_df['igroup']==g+1,'igroup'] = group_name

# Set the first two columns as the MultiIndex
particle_bin_df.set_index(['igroup', 'jbin'], inplace=True)

lines_read += lines_to_read

####################################################
# Display the first few rows of the DataFrame
print(particle_bin_df.head())
print('...')
print(particle_bin_df.tail())

####################################################

lines_to_read = nz

# Read the data into a Pandas DataFrame
atmo_structure_df = pd.read_csv(
    file_path,
    sep='\s+',  # Automatically handles varying spaces as delimiters
    skiprows=lines_read,  # Skip the first line (metadata already read)
    nrows=lines_to_read,  # Read only the required number of lines
    header=None,  # No header in the file
    names=['klevel','z','p','T'],  # Assign column names
    usecols=[0,1,3,4]
)

# Set the first column as the index
atmo_structure_df.set_index(['klevel'], inplace=True)

lines_read += lines_to_read

####################################################
# Display the first few rows of the DataFrame
print(atmo_structure_df.head())
print('...')
print(atmo_structure_df.tail())

####################################################
# Parameters
chunk_size = nz * nbin  # Rows per chunk (metadata + data)

print(len(data_name_list))

data_df = pd.read_csv(
        file_path,
        sep='\s+',  # Automatically handles varying spaces as delimiters
        index_col=False,
        skiprows=lines_read+1,  # Skip the first line (metadata already read)
        #nrows=chunk_size,  # Read only the required number of lines
        header=None,  # No header in the files
        names=data_name_list
)
print(data_df)

print(len(data_df))

####################################################
# Create a boolean mask to keep all rows except every i-th row
mask = [True]*chunk_size
while len(mask) < len(data_df):
    mask+=[False]+[True]*chunk_size

print(len(mask))

notmask = [not m for m in mask]

# Apply the mask
cleaned_df = data_df[mask].copy()
metadata_df = data_df[notmask]

####################################################
#gonna do some janky processing of the metadata_df
timestep = metadata_df['jbin'].values
distance = metadata_df['klevel'].values
rotation = metadata_df['tio2'].values

print(len(timestep))
print(len(distance))
print(len(rotation))

current_longitude = distance - (len(longitude)*rotation)

print(type(current_longitude))
print(len(longitude))

current_longitude = np.round(current_longitude,0)
current_longitude[current_longitude==len(longitude)] = 0

####################################################
#Ok lets add columns
timestep_column_list = [0]*chunk_size
longitude_column_list = [0]*chunk_size
for t,time in enumerate(timestep):
    timestep_column_list.extend([time]*chunk_size)
    longitude_column_list.extend([current_longitude[t]]*chunk_size)

print(len(timestep_column_list))
print(len(longitude_column_list))

cleaned_df['timestep'] = timestep_column_list
cleaned_df['longitude'] = longitude_column_list

#print(cleaned_df)

####################################################
# Set the columns as the index
cleaned_df.set_index(['timestep','jbin','klevel'], inplace=True)

####################################################
# Display the first few rows of the DataFrame
print(cleaned_df.head())
print('...')
print(cleaned_df.tail())


'''# Add columns for the mass and normalising the bins
for g, group_name in enumerate(group_name_list):
    #make a new column
    col_name = 'mass_%s'%group_name
    cleaned_df[col_name] = None
    mass_density = M_dict['M_%s'%group_name]
    for j in range(1,nbin+1):
        print('Doing:',group_name,j)
        thing = cleaned_df.loc[(slice(None), j, slice(None)),group_name]
        bin_radius = particle_bin_df.loc[(group_name,j),'r']
        mass_conv_factor = mass_density*((4./3.)*np.pi*(bin_radius*1e-4)**3.)

        cleaned_df.loc[(slice(None), j, slice(None)),col_name] = thing*mass_conv_factor'''

# Compute mass conversion factors for each group_name
for group_name in group_name_list:
    mass_density = M_dict[f'M_{group_name}']
    particle_bin_df[f'mass_conv_factor'] = (
        particle_bin_df['r'] * 1e-4
    ).apply(lambda r: mass_density * ((4. / 3.) * np.pi * r ** 3.))

    print(particle_bin_df)
quit()

# Reset MultiIndex for merging
merged_df = cleaned_df.reset_index().merge(
    particle_bin_df.reset_index(),  # Bring in particle_bin_df with mass_conv_factor columns
    on=['igroup', 'jbin'],         # Merge on these levels
    how='left'
).set_index(cleaned_df.index.names)  # Restore the original index of cleaned_df

# Multiply and compute mass columns
for group_name in group_name_list:
    cleaned_df[f'mass_{group_name}'] = (
        merged_df[group_name] * merged_df[f'mass_conv_factor_{group_name}']
    )

####################################################
# Display the first few rows of the DataFrame
print(cleaned_df.head())
print('...')
print(cleaned_df.tail())