# Eplink

The current repository stores all codes used for stimulus presentation and analyzing Eplink data (Phase II, Phase III, and TWH) along with a detailed walkthrough of the preprocesing and analysis.

## Datasets

### Phase II

Raw data is available in BIDS format at ComputeCanada:
`/project/6050199/akhanf/cfmm-bids/data/Peters/tle3T_phase2/`

### Phase III

Raw data is available in BIDS format at ComputeCanada:
`/project/6050199/akhanf/cfmm-bids/data/Burneo/EpLinkPhase3_Baseline/`

### TWH

Raw data is available in BIDS format at ComputeCanada:
`project/6050199/akhanf/ext-bids/eplink_phase3`

## Setup

### computecanada

### neuroglia-helpers

### snakemake and snakebids

### wb-command

## Preprocessing

### Step 1: freesurfer 7.2

Datasets were analyzed in **freesurfer** version **7.2**. `regularBatch` from [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) was used to submit recon-all jobs for each participant to compute canada.

**regularBatch usage:** 
```
regularBatch <script_name> <participants_list_text> -b <before-args> -a <after-args> -j <job-template>
```

`regularBatch` takes the following options for passing arguments to pipeline script:
- -b /args/ : cmd-line arguments that go *before* the subject id
- -a /args/ : cmd-line arguments that go *after* the subject id

**For running all subjects:**
```
regularBatch ./run_freesurfer_bids_7.2 subjects.txt -a "\"datasets/eplink-p3/bids eplink-p3-freeSurfer-out/freesurfer_7.2\"" -j 8core32gb12h
```
`run_freesurfer_bids_7.2`: script that runs freesurfer recon-all for one subject (available in this repository).

`subjects.txt`: stores the subjects IDs that is going to be analyzed (included for phase III). You can change this file to rerun freesurfer on arbitrary subjects. One subject id in each line including the 'sub-' (example for all phase III subjects is included in this repository).

`-j 8core32gb12h`: 



### Step 2: fMRIprep 20.2.6
Datasets were preprocessed by the **fMRIprep** pipeline version **20.2.6**. `bidsBatch` from [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) was used to submit fMRIprep jobs for each participant to compute canada. 

**bidsBatch usage:** 
```
bidsBatch <bidsBatch options> <app_name> <bids_dir> <output_dir> <participant/group> <app options>
```

**For running all subjects:**
```
bidsBatch fMRIprep_20.2.6 <bids_dir> <output_dir> participant --output-spaces MNI152NLin2009cAsym T1w --fs-subjects-path <freesurfer_output_dir>
```

`<bids_dir>` and `<output_dir>` can be relative or absolute paths to raw data and output folder. fMRIprep option `--output-spaces` is set to generate outputs in both native (`T1w`) and MNI template spaces (`MNI152NLin2009cAsym`). 

default: OUTPUT_DIR/freesurfer

### Step 3: ISC pipeline

