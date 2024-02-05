[![Static Badge](https://img.shields.io/badge/UWO%20-%20CoNCH?logoColor=%234F2683&label=CoNCH%20lab&labelColor=%238F55E0&color=%234F2683)](https://www.conchlab.uwo.ca)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Eplink

The current repository stores codes used for stimulus presentation and analyzing Eplink data (Phase II, Phase III, and TWH) along with a detailed walkthrough of the preprocesing and analysis.

## Datasets

Raw data is available in BIDS format on Graham (ComputeCanada):

| dataset   | path                                                                   | symbolic links                                     |
|-----------|------------------------------------------------------------------------|----------------------------------------------------|
| Phase II  | `/project/6050199/akhanf/cfmm-bids/data/Peters/tle3T_phase2/`          | `/project/6050199/alit/EpLink/datasets/eplink-p2`  |
| Phase III | `/project/6050199/akhanf/cfmm-bids/data/Burneo/EpLinkPhase3_Baseline/` | `/project/6050199/alit/EpLink/datasets/eplink-p3`  |
| TWH       | `/project/6050199/akhanf/ext-bids/eplink_phase3`                        | `/project/6050199/alit/EpLink/datasets/eplink-twh` |

### Phase II

### Phase III

### TWH


## Setup

### computecanada

### neuroglia-helpers

### snakemake and snakebids

### wb-command

## Preprocessing

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

### Step 3: ISC pipeline (snakemake)

The inter-subject correlation pipline is implemented using the snakemake workflow available in this repository.

#### Setting up python virtual environment

1. Set up a virtual environment: `python3.9 -m venv ~/venv-eplink`
2. Activate the environment: `source ~/venv-eplink/bin/activate`
3. Update pip: `pip install --upgrade pip`
4. Install snakebids: `pip install snakebids`

#### Dry run
Make sure you have activated the virtual environment:

```source ~/venv-eplink/bin/activate```

Change directory to the pipeline folder:

```cd ISC-pipeline```

Dry run does not execute anything, and display what would be done. You can use `--dry-run`, `--dryrun`, or `-n` options to see the missing files and the rules that will make them. You can combine dryrun with `--quiet` or `-q` to just print a summary of the DAG of jobs.

```snakemake -nq```

You can save the extended output of the dry run to a file for examining later:

```snakemake -n > dryrun.txt```

#### Running the pipeline

1. `regularInteractive`
2. Load freesurfer module: `module load freesurfer/7.2.0`
3. `snakemake -c8`

If snakemake stops unexpectedly (e.g., due to job running out of time, or power outage) you'll get the following error indicating that the directory is already locked by an instance of snakemake:

```
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
/project/6050199/alit/EpLink/Eplink/ISC-pipeline
If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
```
As mentioned in the message above, if you are sure that no other snakemake instance is running in this directory you can unlock it by `snakemake --unlock` first and then run snakemake as per usual.
