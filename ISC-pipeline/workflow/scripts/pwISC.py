import nibabel as nib
import numpy as np
import pandas as pd
import re
import h5py
import os
from utils import *

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
# target_volumes = 240
# for original phase III = 384
# for original phase II = 240
# Exclude files with less volumes
df = df[df["#Vols"] >= snakemake.params.n_volumes+snakemake.params.n_dummies].reset_index(drop=True)
# Get unique subjects
subjects = df['Subject'].unique()

# Load files data
data = []
for subj in subjects:
    # Filter subject files 
    df_s = df[df['Subject'] == subj]
    # Load runs
    data.append(load_runs(df_s, snakemake.params.n_volumes, snakemake.params.n_dummies))

# Stacking subjects data (Subject, Unit, Time)
data = np.stack(data)
# Reorganizing data to (Unit, Subject, Time) shape for more efficient slicing when calculating ISCs
data = np.transpose(data, axes=(1, 0, 2))

#### Pair-wise ISC
pw_ISC = [vectorize_lt_matrix(np.corrcoef(data[u,:,:], rowvar=True)) for u in range(data.shape[0])] # we save the lower triangle only
pw_ISC = np.stack(pw_ISC)

# Save ISCs as HDF5 file
with h5py.File(snakemake.output.h5, 'w') as f:
    f.create_dataset('pw_ISC', data=pw_ISC)
    f.create_dataset('subjects', data=subjects)
