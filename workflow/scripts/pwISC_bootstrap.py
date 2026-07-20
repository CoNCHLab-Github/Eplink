import h5py
import numpy as np
from scipy.stats import false_discovery_control as fdr

from utils import recon_lt_matrix


def load_pwISC(fpath):
    """Load a pairwise ISC file, reconstructing matrices from their stored lower triangle."""
    decode = np.vectorize(lambda b: b.decode("utf-8"))

    with h5py.File(fpath, "r") as f:
        pw_ISC_raw = f["pw_ISC"][:]
        if pw_ISC_raw.ndim <= 2:
            pw_ISC = np.stack(
                [recon_lt_matrix(pw_ISC_raw[u]) for u in range(pw_ISC_raw.shape[0])]
            )
        else:
            pw_ISC = pw_ISC_raw
        subjects = decode(f["subjects"][:])
    return pw_ISC, subjects


def bootstrap_pw_matrix(pw_matrix, rng):
    """Shuffle subjects within a similarity matrix, per the recommendation in Chen et al., 2016."""
    n_sub = pw_matrix.shape[0]
    bootstrap_subject = sorted(rng.choice(n_sub, size=n_sub, replace=True))
    return pw_matrix[bootstrap_subject, :][:, bootstrap_subject]


# Concatenate the per-hemisphere pairwise ISC matrices (Unit, Subject, Subject)
ISC = []
subjects = None
for fpath in snakemake.input.func:
    isc, subjects = load_pwISC(fpath)
    ISC.append(isc)
pw_ISC = np.concatenate(ISC, axis=0)

# Restrict to control subjects, excluding any flagged for QC
subjects = np.array(list(map(int, subjects)))
controls = subjects > 5000
controls[np.isin(subjects, snakemake.params.subjid_exclude)] = False
pw_ISC = pw_ISC[:, controls, :][:, :, controls]
pw_ISC[np.isnan(pw_ISC)] = 0

n_unit, n_sub, _ = pw_ISC.shape
n_perm = snakemake.params.n_perm
rng = np.random.default_rng(snakemake.params.seed)

# Subject-Wise Bootstrapping (SWB): resample the similarity matrix n_perm times per unit and
# test whether the median ISC is significantly greater than the null (median < 0)
# PERF TODO: this is O(n_unit * n_perm) individual np.median calls, which is prohibitively slow
# for vertex-level (atlas=none) units. A shared-resample-per-permutation approach (draw one
# bootstrap subject index array per permutation, gather+median across ALL units at once with
# pw_ISC[:, idx, :][:, :, idx]) cuts this to O(n_perm) numpy calls and should be revisited once
# we can benchmark properly on a compute node (not the login node).
stat = np.zeros((n_unit, n_perm))
for u in range(n_unit):
    for p in range(n_perm):
        stat[u, p] = np.median(bootstrap_pw_matrix(pw_ISC[u], rng))

pval = (stat < 0).mean(axis=1)
pval_adj = fdr(pval, method="by")

with h5py.File(snakemake.output.h5, "w") as f:
    f.create_dataset("pval", data=pval)
    f.create_dataset("pval_adj", data=pval_adj)
    f.create_dataset("controls", data=controls)
