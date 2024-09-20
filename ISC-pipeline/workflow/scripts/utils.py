import nibabel as nib
import numpy as np
import pandas as pd
import re
import h5py
import os

def get_dimensions(filepath):
    """Extract dimensions of data from h5 data."""
    # Check file extension
    ext = os.path.splitext(filepath)[-1]
    
    # h5 file
    if ext == '.h5': 
        # Open the HDF5 file
        with h5py.File(filepath, 'r') as f:
            # Load the parcellated data
            data = f['parcellated_data'][:]
            # Return data dimensions
            n_rois, n_vol = data.shape
    
    # gii file
    elif ext == '.gii':
        # Load the GIFTI file by nibabel
        gii = nib.load(filepath)
        # Return data dimensions
        n_rois = gii.darrays[0].data.shape[0]
        n_vol = len(gii.darrays)
    
    return n_rois, n_vol

def load_HDF(filepath, n_vols, n_dummies=0):
    """Load data stored in the HDF files given the runs dataframe."""
    with h5py.File(filepath, 'r') as f:
            # Load the parcellated data
            data = f['parcellated_data'][:]
            data = data[:,n_dummies:n_dummies+n_vols] # Trimming excessive volumes
    # Returns data with shape ROI x Time
    return data

def load_gii(filepath, n_vols, n_dummies=0):
    func_gii = nib.load(filepath)
    data = np.vstack([darray.data for darray in func_gii.darrays[n_dummies:n_dummies+n_vols]]).T
    # Returns data with shape Vertex x Time
    return data

def load_runs(runs_df, n_vols, n_dummies=0):
    # Loading only the last run
    fp = runs_df['Path'].iloc[-1]
    # Checking file extension
    ext = os.path.splitext(fp)[-1]

    # h5 file
    if ext == '.h5':
        data = load_HDF(fp, n_vols, n_dummies)
    # gii file
    elif ext == '.gii':
        data = load_gii(fp, n_vols, n_dummies)
    
    return data

def vectorize_lt_matrix(pw_matrix:np.ndarray):
    return pw_matrix[np.tril_indices_from(pw_matrix, k=-1)]

def recon_lt_matrix(pw_vec:np.ndarray, n=None):
    if n is None:
        n = int(1+np.sqrt(1+8*pw_vec.size))//2
    pw_matrix = np.zeros((n,n))
    pw_matrix[np.tril_indices(n,k=-1)] = pw_vec
    pw_matrix = pw_matrix+pw_matrix.T
    pw_matrix[np.diag_indices(n)] = 1
    return pw_matrix