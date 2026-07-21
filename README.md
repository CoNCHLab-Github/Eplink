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
    - [Step 3: ISC pipeline (snakemake + SLURM, via pixi)](#step3)

## Datasets

Raw data and derivatives are available in BIDS format on Nibi (Alliance Canada), under the CoNCH lab's allocation:

| dataset            | path                                        | notes                                                                                                    |
|--------------------|----------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| EpLink3 (Phase III) | `/project/6033493/datasets/EpLink3`         | `dataset_description.json` calls it `EpLink_merged`: combines the original EpLink3 cohort (syngo MR E11) with an additional healthy-control cohort, EpLink_newHC (syngo MR XA60). Cohort is recorded per-subject in `participants.tsv`. freesurfer/fmriprep derivatives already present as per-subject tar.gz (see below). |
| EpLink2 (Phase II)  | `/project/6033493/datasets/EpLink2/rawdata` | Raw BIDS data only; freesurfer/fmriprep derivatives have not yet been generated under this allocation.  |

Each dataset directory follows the same layout, which the pipeline's config templates against (`config/snakebids.yml`):

```
/project/6033493/datasets/<dataset>/
├── participants.tsv
├── dataset_description.json
└── derivatives/
    ├── freesurfer-7.2.0/sub-<id>_freesurfer-7.2.0.tar.gz
    └── fmriprep-20.2.6/sub-<id>_fmriprep-20.2.6.tar.gz
```

Which dataset is processed is controlled entirely by the `dataset:` key in `config/snakebids.yml` (currently `EpLink3`) — it gets substituted into every `{dataset}`-templated path in the config (`participants_tsv`, `freesurfer_tar`, `fmriprep_tar`, `archive_dir`, ...).

### Phase II

This data is collected from epileptic patients, and age and sex matched healthy controls at the Centre for Functional and Metabolic Mapping, Western University, London, Ontario. The dataset consists of a structural scan and three naturalistic tasks: a resting state and a movie watching task. The movie was *Bang! You're Dead* by Alfred Hitchcock (hitchcock). The following table summarizes the acquisition parameters for the tasks:

| run       | # Volumes | TR (ms)   | TE (ms) | Slice Thickness (mm) |
| --------- | --------- | --------- | ------- | -------------------- |
| rest      | 175       | 2200      | 30      | 3.7                  |
| hitchcock | 246       | 2000      | 30      | 3                    |

### Phase III

### TWH

## Setup

### 1. Get access and clone the repository

The pipeline reads datasets from the CoNCH lab's Alliance Canada allocation (`/project/6033493/...`, account `rrg-conchlab_cpu`) and submits SLURM jobs on Nibi — you need an active Alliance Canada account with access to that allocation before anything else here will work (ask the PI if you don't have this yet).

Clone the repo on Nibi (e.g. into your `/scratch/<username>` space, since that's where the pipeline's `results/` outputs land):
```
git clone https://github.com/CoNCHLab-Github/Eplink.git
```

### 2. computecanada

#### checking the submitted jobs and their status

use `squeue` or `sq` command to list Slurm jobs. `squeue` supplies information about all jobs in the system, by default. You can use the shorter `sq` to list only your own jobs. Find more details [here](https://docs.alliancecan.ca/wiki/Running_jobs#Use_squeue_or_sq_to_list_jobs).

<div id="sshfs"/>

#### mounting compute canada as a file system using sshfs:
1. create a directory for the mount if you don't have one: `mkdir ~/Nibi`
2. use sshfs to mount as a directory on computecanada file system: ```sshfs <username>@nibi.sharcnet.ca:/home/<username> ~/Nibi```
3. You will be asked for your password and 2FA if applicable.

### 3. pixi

The ISC pipeline's Python environment (snakemake, snakebids, nilearn, etc.) and its run commands are managed by [pixi](https://pixi.sh) — dependencies are pinned in `pixi.toml`/`pixi.lock` so the environment is reproducible across machines.

1. Install pixi once (works on Nibi or locally): `curl -fsSL https://pixi.sh/install.sh | bash`
2. From the repo root, create the environment: `pixi install`
3. Run any pipeline command through pixi — it resolves and activates the environment automatically, no manual `source activate`/`venv` step needed:
   ```
   pixi run <task>
   ```
4. See the full list of available tasks with `pixi task list`, or read the `[tasks]` table in `pixi.toml`. Tasks are described in [Step 3](#step3) below.
5. An optional `dev` environment (`pytest`, `ipython`, `jupyterlab`) is available via `pixi run -e dev <task>`.

### 4. wb-command

`wb_command` (Connectome Workbench) is not installed via pixi. On Nibi it's loaded as an environment module — the module name is set by `wb_module` in `config/snakebids.yml` (default: `connectome-workbench`; find the exact name with `module spider connectome`). Rules that call `wb_command` also load `apptainer` alongside it (some `wb_command` operations shell out to a containerized helper), and Snakemake loads both modules automatically for SLURM-submitted jobs because the profile sets `use-envmodules: true` (see [Step 3](#step3)).

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

<details>
<summary>

**Example command for phase II**
</summary>

```
regularBatch ./run_freesurfer_bids_7.2 subjects.txt -a "\"/project/6050199/akhanf/cfmm-bids/data/Peters/tle3T_phase2/bids /scratch/alit/eplink-p2-freesurfer\"" -j 8core32gb12h
```
</details>

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
To pull the fMRIprep reports properly on your local machine to review run the `fetch_fMRIprep_reports.sh` script (included in the repository). You might need to adjust the fMRIprep output path in the script. This script assumes that Nibi is mounted as a file system at `~/Nibi`, if it isn't check [here](#sshfs) for instructions.

---
<div id="step3"/>

### Step 3: ISC pipeline (snakemake + SLURM, via pixi)

The inter-subject correlation pipeline is implemented as a [snakebids](https://snakebids.readthedocs.io)-wrapped Snakemake workflow: `workflow/Snakefile` (rules) + `config/snakebids.yml` (dataset selection, parcellations, confounds, smoothing levels, etc.). All commands are run through `pixi run <task>` (see [pixi setup](#setup)), which activates the pinned environment for you — no manual venv/module activation needed for the Python side.

```
pixi install                     # once, creates the pinned environment
pixi run dry-run                 # sanity-check the DAG (snakemake -n)
pixi run extract-fmriprep        # one-time: untar fmriprep derivatives (see below)
pixi run run                     # submit the rest of the pipeline to SLURM
```

#### How jobs get submitted to SLURM

`workflow/profiles/slurm/config.yaml` is a Snakemake 8 [executor profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) built on `snakemake-executor-plugin-slurm` (pinned in `pixi.toml`). Every `pixi run` pipeline task passes `--profile workflow/profiles/slurm` to snakemake, so when you invoke it snakemake:

- resolves the DAG of missing output files for the target rule,
- submits each job individually to SLURM (`sbatch`) under `slurm_account: rrg-conchlab_cpu`,
- applies default resources (`mem_mb: 32000`, `cpus_per_task: 8`, `runtime: 60` min) unless a rule overrides them in its own `resources:` block in the Snakefile,
- loads environment modules (e.g. `wb_command`) for rules that declare `envmodules:`, since the profile sets `use-envmodules: true`,
- retries failed jobs twice (`retries: 2`) and keeps unrelated jobs running on failure (`keep-going: true`),
- runs up to 100 jobs concurrently (`jobs: 100`).

You do not write or submit an `sbatch` script yourself — run `pixi run run` (or any of the stage-specific tasks below) from a login node or inside `regularInteractive`, and Snakemake handles submission, polling, and resubmission for every rule. Use `squeue`/`sq` (see [above](#computecanada)) to watch the submitted jobs.

#### Pipeline stages and matching pixi tasks

Each stage has a `dry-<stage>` (adds `--dry-run`) and a `<stage>` (submits for real) pixi task:

| stage | snakemake target | pixi tasks | what it does | output path |
|---|---|---|---|---|
| fMRIprep extraction | `extract_fmriprep_all` | `dry-extract-fmriprep`, `extract-fmriprep` | untars each subject's `fmriprep-20.2.6` derivatives tar into `results/{dataset}/fmriprep/`. **Must run once before the rest of the pipeline** — later rules read from this extracted directory. | `results/{dataset}/fmriprep/` |
| denoising | `denoise_all` | `dry-denoise`, `denoise` | nilearn-based confound regression on the T1w-space preprocessed BOLD (`workflow/scripts/denoise.py`), per confound set in `config["confounds"]`. | `results/{dataset}/denoised/` |
| surfaces | `gen_midthickness_all` | `dry-surfaces`, `surfaces` | extracts freesurfer surfaces from `freesurfer_tar`, converts white/pial to GIFTI directly in scanner-RAS space (via `mris_convert --to-scanner`), and builds native-space midthickness surfaces (needs the freesurfer container + `wb_command`). | `results/{dataset}/surfaces/` |
| resampling | `resample2fsLR_all` | `dry-fslr`, `fslr` | projects denoised volumes onto the native surface, then resamples to fsLR 32k. | `results/{dataset}/resampled2fsLR32k/` |
| smoothing | `smooth_all` | `dry-smooth`, `smooth` | surface smoothing at each `fwhm` in `config["fwhm"]` (FWHM is converted to sigma before being passed to `wb_command`). | `results/{dataset}/smoothed/` |
| parcellation | `parcellate_original_all` | `dry-parcellate`, `parcellate` | averages resampled (or temporally-resampled) surface data within each `config["parcellations"]` atlas. | `results/{dataset}/{source}-parcellated/{atlas}/` |
| ISC | `ISC_all` | `dry-isc`, `isc` | pairwise and leave-one-out ISC (vertex-level and per-parcellation) plus functional connectivity. | `results/{dataset}/pwISC/`, `results/{dataset}/FC/{atlas}/` |
| ISC significance (SWB) | `pwISC_bootstrap_movie` | `dry-pwisc-swb`, `pwisc-swb` | subject-wise bootstrap significance testing (Chen et al., 2016) on the pairwise ISC from `pw_ISC`, run over every parcellation atlas and every smoothing level (vertex-level included), excluding QC-flagged subjects (`config["subjid_exclude"]`). Vertex-level (`atlas=none`) jobs get a longer SLURM runtime since the bootstrap loop isn't yet optimized for that many units — see the `PERF TODO` in `workflow/scripts/pwISC_bootstrap.py`. | `results/{dataset}/pwISC_stats/` |
| archiving | `archive_results` | `dry-archive`, `archive` | packages `results/{dataset}` into a SquashFS image under `archive_dir` (the dataset's project-space directory). | `archive_dir` (project space) |
| everything | (default target) | `dry-run`, `run` | runs the full DAG up to the pipeline's default target. | — |
| DAG visualization | — | `dag` | `snakemake --dag \| dot -Tsvg > dag.svg` | — |

#### Configuring the pipeline

**Parcellations**

Here is a list of available parcellations:
- [Desikan](https://doi.org/10.1016/j.neuroimage.2006.01.021)
- [Glasser 2016](https://doi.org/10.1038%2Fnature18933)
- [Schaefer 2018](https://doi.org/10.1093/cercor/bhx179)
- [Yan 2023](https://doi.org/10.1016/j.neuroimage.2023.120010)

> **Note 1:** Desikan, Glasser, and Schaefer atlas files were adopted form [this repository](https://github.com/DiedrichsenLab/fs_LR_32).

> **Note 2:** Yan 2023 parcellations come in CIFTI format and were converted to .gii using wb_command:
`
wb_command -cifti-separate <input>.dlabel.nii COLUMN -label CORTEX_LEFT <output_L>.label.gii -label CORTEX_RIGHT <output_R>.label.gii
`

**ISC significance testing (SWB)**

`n_perm_bootstrap` (`config/snakebids.yml`) sets the number of subject-wise bootstrap resamples used by `pwISC_bootstrap` (default 10000); `bootstrap_seed` fixes the RNG seed for reproducibility (leave blank for a random seed per run). `subjid_exclude` lists subject IDs to drop from the significance test (e.g. QC failures), matching the exclusions applied in `notebooks/pairwise_summary.ipynb`.

#### Troubleshooting: locked directory

If snakemake stopped unexpectedly (e.g., a job ran out of time, or a power outage) you'll get an error indicating the run directory is locked by a previous instance:

```
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
...
If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
```

If you're sure no other snakemake instance is running, unlock it first, then rerun as usual:

```
pixi run snakemake --snakefile workflow/Snakefile --profile workflow/profiles/slurm --unlock
```
