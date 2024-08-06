import numpy as np
import h5py
import os

def load_HDF(filepath, n_vols=None):
    """Load data stored in the HDF files given the runs dataframe."""
    if n_vols == None:
        n_vols = -1
    with h5py.File(filepath, 'r') as f:
            # Load the parcellated data
            data = f['parcellated_data'][:]
            data = data[:,:n_vols] # Ignoring excessive volumes
    # Returns data with shape ROI x Time
    return data

data = []
for f in snakemake.input.func:
    data.append(load_HDF(f))

# Data with shape ROIs x Time
data = np.vstack(data)

pearson_FC = np.corrcoef(data, rowvar=True)

# Save FC matrix as HDF5 file
with h5py.File(snakemake.output.h5, 'w') as f:
    f.create_dataset('pearson_FC', data=pearson_FC)
