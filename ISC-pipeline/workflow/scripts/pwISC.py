import nibabel as nib
import numpy as np
import h5py

# f = open(snakemake.output.h5, "a")
# f.write(type(snakemake.input.func))
# f.close()
for f in snakemake.input.func:
    print(f)
    #f is string

# # Load functional data and atals file
# func_gii = nib.load(snakemake.input.func)
# label_gii = nib.load(snakemake.input.atlas)

# # Extract the data arrays from the loaded GIfTI files
# func_data = np.vstack([darray.data for darray in func_gii.darrays])
# label_data = label_gii.darrays[0].data
# label_table = label_gii.labeltable.get_labels_as_dict()

# # Get parcels from the label table
# ROIs = label_table.keys()

# # Initialize a list to hold the parcellated functional data
# parcellated_data = np.zeros((len(ROIs),func_data.shape[1]))

# # Parcellate the functional data
# for roi in ROIs:
#     # Find vertices that belong to the current parcel
#     vertices_in_parcel = np.where(label_data == roi)[0]
    
#     # Extract the functional data for these vertices and average across vertices
#     parcel_time_series = np.mean(func_data[vertices_in_parcel, :], axis=0)
    
#     # Append the averaged time series to the parcellated data list
#     parcellated_data[roi,:] = parcel_time_series

# # Save Parcelated data as HDF5 file
# with h5py.File(snakemake.output.h5, 'w') as f:
#     f.create_dataset('parcellated_data', data=parcellated_data)