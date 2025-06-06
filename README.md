[![Static Badge](https://img.shields.io/badge/UWO%20-%20CoNCH?logoColor=%234F2683&label=CoNCH%20lab&labelColor=%238F55E0&color=%234F2683)](https://www.conchlab.uwo.ca)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Eplink

The current repository stores codes used for stimulus presentation and analyzing Eplink data (Phase II, Phase III, and TWH) along with a detailed walkthrough of the preprocesing and analysis.

**Table of Contents**
- [Datasets](#datasets)
- [Setup](#setup)
- [Preprocessing](#preprocessing)
    - [Step 1: freesurfer 7.2](#step1)
    - [Step 2: fMRIprep 20.2.6](#step2)
    - [Step 3: ISC pipeline (snakemake)](#step3)

## Datasets

Raw data is available in BIDS format on Graham (ComputeCanada):

| dataset   | path                                                                       | symbolic links                                     |
|-----------|----------------------------------------------------------------------------|----------------------------------------------------|
| Phase II  | `/project/6050199/khanlab/datasets/internal/Peters_EpLink/tle3T_phase2/`   | `/project/6050199/alit/EpLink/datasets/eplink-p2`  |
| Phase III | `/project/6050199/khanlab/datasets/internal/Burneo_EpLinkPhase3_Baseline/` | `/project/6050199/alit/EpLink/datasets/eplink-p3`  |
| TWH       | `/project/6050199/khanlab/datasets/external/eplink_phase3`                 | `/project/6050199/alit/EpLink/datasets/eplink-twh` |

### Phase II

This data is collected from epileptic patients, and age and sex matched healthy controls at the Centre for Functional and Metabolic Mapping, Western University, London, Ontario. The dataset consists of a structural scan and three naturalistic tasks: a resting state and a movie watching task. The movie was *Bang! You're Dead* by Alfred Hitchcock (hitchcock). The following table summarizes the acquisition parameters for the tasks:

| run       | # Volumes | TR (ms)   | TE (ms) | Slice Thickness (mm) |
| --------- | --------- | --------- | ------- | -------------------- |
| rest      | 175       | 2200      | 30      | 3.7                  |
| hitchcock | 246       | 2000      | 30      | 3                    |

### Phase III

### TWH

## Setup

### computecanada

#### checking the submitted jobs and their status

use `squeue` or `sq` command to list Slurm jobs. `squeue` supplies information about all jobs in the system, by default. You can use the shorter `sq` to list only your own jobs. Find more details [here](https://docs.alliancecan.ca/wiki/Running_jobs#Use_squeue_or_sq_to_list_jobs).

<div id="sshfs"/>

#### mounting compute canada as a file system using sshfs:
1. create a directory for the mount if you don't have one: `mkdir ~/Graham`
2. use sshfs to mount as a directory on computecanada file system: ```sshfs <username>@graham.computecanada.ca:/home/<username>/projects/ctb-akhanf/<username> ~/Graham```
3. You will be asked for your password and 2FA if applicable.

### neuroglia-helpers

### snakemake and snakebids

### wb-command

## Preprocessing

<div id="step1"/>

### Step 1: freesurfer 7.2 

Datasets were analyzed in **freesurfer** version **7.2**. `regularBatch` from [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) was used to submit recon-all jobs for each participant to compute canada.

#### regularBatch usage:
```
regularBatch <script_name> <participants_list_text> -b <before-args> -a <after-args> -j <job-template>
```

`regularBatch` takes the following options for passing arguments to pipeline script:
- -b /args/ : cmd-line arguments that go *before* the subject id
- -a /args/ : cmd-line arguments that go *after* the subject id

#### For running subjects in `subjects.txt`:
```
regularBatch ./run_freesurfer_bids_7.2 subjects.txt -a "\"path/to/bids path/to/output\"" -j 8core32gb12h
```
- `run_freesurfer_bids_7.2`: script that runs freesurfer recon-all for one subject (available in this repository).
- `subjects.txt`: stores the subjects IDs that is going to be analyzed. You can change this file to rerun freesurfer on arbitrary subjects. One subject id in each line including the 'sub-' (example for all phase III subjects is included in this repository).
- `-a "\"path/to/bids path/to/output\""`: passes both bids and output directories to the script as one argument (space seperated).
- `-j 8core32gb12h`: specifies job template for each participant (8 cores, 32GB RAM, 12 hours)


---
<div id="step2"/>

### Step 2: fMRIprep 20.2.6
Datasets were preprocessed by the **fMRIprep** pipeline version **20.2.6**. `bidsBatch` from [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) was used to submit fMRIprep jobs for each participant to compute canada. 

#### bidsBatch usage:
```
bidsBatch <bidsBatch options> <app_name> <bids_dir> <output_dir> <participant/group> <app options>
```
- For the full list of available apps run `bidsBatch` command.
- Find a full list of fMRIprep options [here](https://fmriprep.org/en/stable/usage.html).

#### For running all subjects:
```
bidsBatch fmriprep_20.2.6 <bids_dir> <output_dir> participant --output-spaces MNI152NLin2009cAsym T1w --fs-subjects-dir <freesurfer_output_dir>
```

- `<bids_dir>` and `<output_dir>` can be relative or absolute paths to raw data and output folder.  
- fMRIprep option `--output-spaces` is set to generate outputs in both native (`T1w`) and MNI template spaces (`MNI152NLin2009cAsym`).
- `--fs-subjects-dir` is the path to freesurfer output (default: OUTPUT_DIR/freesurfer).

#### Reviewing fMRIprep reports:
To pull the fMRIprep reports properly on your local machine to review run the `fetch_fMRIprep_reports.sh` script (included in the repository). You might need to adjust the fMRIprep output path in the script. This script assumes that graham is mounted as a file system at `~/Graham`, if it isn't check [here](#sshfs) for instructions.

---
<div id="step3"/>

### Step 3: ISC pipeline (snakemake)

The inter-subject correlation pipline is implemented using the snakemake workflow available in this repository.

![Alt text](./ISC-pipeline/dag.svg)

#### Setting up python virtual environment

1. Set up a virtual environment: `python3.9 -m venv ~/venv-eplink`
2. Activate the environment: `source ~/venv-eplink/bin/activate`
3. Update pip: `pip install --upgrade pip`
4. Install dependencies: `pip install -r requirements.txt`

You can find the `requirements.txt` file in this repository. 

#### Dry run
Make sure you have activated the virtual environment:

```source ~/venv-eplink/bin/activate```

Change directory to the pipeline folder:

```cd ISC-pipeline```

Dry run does not execute anything, and display what would be done. You can use `--dry-run`, `--dryrun`, or `-n` options to see the missing files and the rules that will make them. You can combine dryrun with `--quiet` or `-q` to just print a summary of the DAG of jobs.

```snakemake -nq```

You can save the extended output of the dry run to a file for examining later:

```snakemake -n > dryrun.txt```

#### Generate a dag of rules (visualize the pipeline) 
```snakemake --forceall --rulegraph | dot -Tsvg > dag.svg```

#### Running the pipeline

1. `regularInteractive`
2. Load freesurfer module: `module load freesurfer/7.2.0`
3. `snakemake -c8`

If snakemake had stopped unexpectedly (e.g., due to job running out of time, or power outage) you'll get the following error indicating that the directory is already locked by an instance of snakemake:

```
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
/project/6050199/alit/EpLink/Eplink/ISC-pipeline
If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
```
As mentioned in the message above, if you are sure that no other snakemake instance is running in this directory you can unlock it by `snakemake --unlock` first, and then run snakemake as per usual.
