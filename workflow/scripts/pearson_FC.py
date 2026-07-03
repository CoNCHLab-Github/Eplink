import numpy as np
import h5py
import os
from utils import load_HDF, vectorize_lt_matrix

data = []
for f in snakemake.input.func:
    data.append(load_HDF(f, snakemake.params.n_volumes, snakemake.params.n_dummies))

# Data with shape ROIs x Time
data = np.vstack(data)

pearson_FC = np.corrcoef(data, rowvar=True)
FC_vec = vectorize_lt_matrix(pearson_FC)
params = {'n_vols':snakemake.params.n_volumes, 'n_dummies':snakemake.params.n_dummies}

# Save FC matrix as HDF5 file
with h5py.File(snakemake.output.h5, 'w') as f:
    f.create_dataset('pearson_FC', data=FC_vec)
