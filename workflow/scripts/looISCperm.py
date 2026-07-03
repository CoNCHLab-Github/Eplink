import nibabel as nib
import numpy as np
import pandas as pd
import re
import h5py
import os

# Helper functions
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

def load_HDF(filepath, n_vols):
    """Load data stored in the HDF files given the runs dataframe."""
    with h5py.File(filepath, 'r') as f:
            # Load the parcellated data
            data = f['parcellated_data'][:]
            data = data[:,:n_vols] # Ignoring excessive volumes
    # Returns data with shape ROI x Time
    return data

def load_gii(filepath, n_vols):
    func_gii = nib.load(filepath)
    data = np.vstack([darray.data for darray in func_gii.darrays[:n_vols]]).T
    # Returns data with shape Vertex x Time
    return data

def load_runs(runs_df, n_vols):
    # Loading only the last run
    fp = runs_df['Path'].iloc[-1]
    # Checking file extension
    ext = os.path.splitext(fp)[-1]

    # h5 file
    if ext == '.h5':
        data = load_HDF(fp, n_vols)
    # gii file
    elif ext == '.gii':
        data = load_gii(fp, n_vols)
    
    return data

# Infer task form the output name
pattern = re.compile(r"task-([^_]+)")
match = pattern.search(snakemake.output.h5)
target_task, = match.groups()

# Gather info about the inputs and build a dataframe
pattern = re.compile(r"sub-(\d+)/.*task-([^_]+)_run-(\d+)")

inputs = []
for f in snakemake.input.func:
    match = pattern.search(f)
    if match:
        subject, task, run = match.groups()
        if task == target_task: # filter for the task of interest
            n_unit, n_vol = get_dimensions(f)
            inputs.append({"Subject": subject, "Task": task, "Run": run, "#ROIs": n_unit, "#Vols": n_vol, "Path": f})

# Build a dataframe of input files
df = pd.DataFrame(inputs)
# Sort files 
df = df.sort_values(by=['Subject', 'Task', 'Run'], ascending=[True]*3).reset_index()
# Target volumes ===> TODO: a method to automatically determining the number of volumes
target_volumes = 240# 384
# Exclude files with less volumes
df = df[df["#Vols"] >= target_volumes].reset_index(drop=True)
# Get unique subjects
subjects = df['Subject'].unique()
n_subject = len(subjects)

# Load files data
data = []
for subj in subjects:
    # Filter subject files 
    df_s = df[df['Subject'] == subj]
    # Load runs
    data.append(load_runs(df_s, target_volumes))

# Stacking subjects data (Subject, Unit, Time)
data = np.stack(data)
# Reorganizing data to (Unit, Subject, Time) shape for more efficient slicing when calculating ISCs
# data = np.transpose(data, axes=(1, 0, 2))

# Sum subjects data (Unit, Time)
data_sum = data.sum(axis=0)

# Output will be shaped as (Subject, Unit)
loo_ISC = np.zeros((n_subject, n_unit))
for i in range(n_subject):
    for j in range(n_unit):
        sum_ts = data_sum[j,:].reshape((1,-1))
        sub_ts = data[i,j,:].reshape((1,-1))
        c = np.corrcoef(sub_ts, sum_ts-sub_ts)
        loo_ISC[i,j] = c[0,1]

# Save ISCs as HDF5 file
with h5py.File(snakemake.output.h5, 'w') as f:
    f.create_dataset('loo_ISC', data=loo_ISC)
