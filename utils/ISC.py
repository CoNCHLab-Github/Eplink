import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

def calculate_ISC(data, subjects:range, method='pairwise-mean', n_g1=None, verbose=True):
    # Calculate ISC
    n_subjects = len(subjects)
    mask = np.tril(np.ones(n_subjects, dtype=bool), k=-1)

    # Initialization
    n_rois = data['L'].shape[0]
    if ((method == 'pairwise-mean') |
        (method == 'pairwise-mat') |
        (method == 'loo')):
        ISC = {'L':[], 'R':[]}
    elif method == 'pairwise-mixed':
        ISC = {'b_L':[], 'b_R':[], # contain between groups correlattions
               'w_g1_L':[], 'w_g1_R':[], # contain within group1 correlattions
               'w_g2_L':[], 'w_g2_R':[] # contain within group2 correlattions
               }
        mask_b = np.zeros((len(subjects),)*2, dtype=bool)
        mask_b[n_g1:, :n_g1] = 1
        mask_wg1 = np.zeros((len(subjects),)*2, dtype=bool)
        mask_wg1[:n_g1,:n_g1] = 1
        mask_wg1 = mask_wg1 & mask
        mask_wg2 = np.zeros((len(subjects),)*2, dtype=bool)
        mask_wg2[n_g1:,n_g1:] = 1
        mask_wg2 = mask_wg2 & mask
    elif ((method == 'Xloo') |
          (method == 'Xloo2')):
        ISC = {'CL':[], 'CR':[],
               'XL':[], 'XR':[]
               }
        
    for h in ['L','R']:
        # verbose
        if verbose:
            print(f'Processing hemisphere: {h}')
        for roi in tqdm(range(n_rois), disable=not verbose):
            corr_mat = np.corrcoef(data[h][roi,subjects,:])
            if method == 'pairwise-mean':
                ISC[h].append(corr_mat[mask].mean())
            elif method == 'pairwise-mat':
                ISC[h].append(corr_mat)
            elif method == 'pairwise-mixed':
                ISC[f'b_{h}'].append(corr_mat[n_g1:, :n_g1]) # dimensions: n_g2 x n_g1
                ISC[f'w_g1_{h}'].append(corr_mat[mask_wg1])
                ISC[f'w_g2_{h}'].append(corr_mat[mask_wg2])
            elif method == 'loo':
                for s in subjects:
                    subjects_loo = [x for x in subjects if x != s]
                    subject_int = data[h][roi,s,:]
                    subjects_mean = np.mean(data[h][roi,subjects_loo,:], axis=0)
                    corr = np.corrcoef(subject_int, subjects_mean)
                    ISC[h].append(corr)
            elif method == 'Xloo':
                g2_mean = np.mean(data[h][roi,range(n_g1,n_subjects),:], axis=0)
                isc = []
                for s in range(n_g1):
                    subject_int = data[h][roi,s,:]
                    corr = np.corrcoef(subject_int, g2_mean)
                    isc.append(corr[0,1])
                ISC[f'X{h}'].append(isc)
                isc = []
                for s in range(n_g1,n_subjects):
                    subjects_loo = [x for x in range(n_g1,n_subjects) if x != s]
                    subject_int = data[h][roi,s,:]
                    subjects_mean = np.mean(data[h][roi,subjects_loo,:], axis=0)
                    corr = np.corrcoef(subject_int, subjects_mean)
                    isc.append(corr[0,1])
                ISC[f'C{h}'].append(isc)
            elif method == 'Xloo2':
                isc = []
                for s in range(n_g1): # iterate on subjects in group 1
                    temp = []
                    subject_g1 = data[h][roi,s,:]
                    for s2 in range(n_g1,n_subjects):
                        # leave one subject out from group 2
                        subjects_loo = [x for x in range(n_g1,n_subjects) if x != s2]
                        g2_subjects_mean = np.mean(data[h][roi,subjects_loo,:], axis=0)
                        
                        corr = np.corrcoef(subject_g1, g2_subjects_mean)
                        temp.append(corr[0,1])
                    isc.append(temp)
                ISC[f'X{h}'].append(isc)
                isc = []
                for s in range(n_g1,n_subjects):
                    subjects_loo = [x for x in range(n_g1,n_subjects) if x != s]
                    subject_int = data[h][roi,s,:]
                    subjects_mean = np.mean(data[h][roi,subjects_loo,:], axis=0)
                    corr = np.corrcoef(subject_int, subjects_mean)
                    isc.append(corr[0,1])
                ISC[f'C{h}'].append(isc)
                    

    return {key: np.array(value) for key, value in ISC.items()}

def _nonparam_test(data, n_perm):
    # Data with dimensions: Subject x Time
    n_subj, n_tp = data.shape
    # Declare return matrix with size Subject x Permutation+1
    isc = np.zeros((n_subj, n_perm+1))
    for lo_sub in range(n_subj):
        # Subject indices excluding the left out subject 
        row_indices = np.array([s for s in range(n_subj) if s!=lo_sub])[:, np.newaxis]
        # Define variable to store original and shifted samples
        shifted_perms = np.zeros((n_perm+1, n_tp))
        # Original synchrony
        shifted_perms[0,:] = np.mean(data[row_indices, :], axis=0)
        # Produce shifted samples 
        for i in range(n_perm):
            shift_values = np.random.randint(0, n_tp, size=n_subj-1)
            col_indices = (np.arange(n_tp) - shift_values[:, np.newaxis])
            shifted_perms[i+1,:] = np.mean(data[row_indices, col_indices], axis=0)
        # Calculate and store the correlations
        corr = np.corrcoef(np.vstack((data[lo_sub, :],shifted_perms)))
        isc[lo_sub,:] = corr[1:,0]
        
    return isc


def nonparam_test(data, n_perm):
    result = dict()
    for hemi in ['L', 'R']:
        print(f'processing {hemi} hemisphere')
        # Data with dimensions: Units x Subjects x Time
        n_units, n_subj, _ = data[hemi].shape

        partial_nonparam_test = partial(_nonparam_test, n_perm=n_perm)
        units = range(n_units)

        result[hemi]  = np.zeros((n_units,n_subj,n_perm+1))

        with ThreadPoolExecutor(max_workers=16) as executor:
            # Map the function to the items and keep track of the futures
            futures = {executor.submit(partial_nonparam_test, data[hemi][unit]): unit for unit in units}

            # Progress bar setup
            with tqdm(total=len(units)) as progress_bar:
                for future in as_completed(futures):
                    # Result of a completed task
                    unit_id = futures[future]  # Get the unit identifier
                    result[hemi][unit_id,:,:] = future.result()
                    # Update the progress bar
                    progress_bar.update(1)

    return result