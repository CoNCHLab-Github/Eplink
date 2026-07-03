import nibabel as nib
import numpy as np
import os
import h5py
from utils import load_HDF

fp = snakemake.input.ts
ext = os.path.splitext(fp)[-1]

if ext == '.h5':
    with h5py.File(fp, 'r') as f:
        # Load the parcellated data
        data = f['parcellated_data'][:] # ROI x Time
    
    tSNR = np.divide(data.mean(axis=1),data.std(axis=1))

    with h5py.File(snakemake.output.tSNR, 'w') as f:
        # Save the tSNR data
        f.create_dataset('tSNR', data=tSNR)

elif ext == '.gii':
    print('not impelented yet')
elif ext == '.nii':
    print('not impelented yet')
else:
    print('Invalid extention')