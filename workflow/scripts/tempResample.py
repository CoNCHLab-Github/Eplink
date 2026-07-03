import nibabel as nib
import numpy as np
#from snakebids import write_derivative_json
from scipy.signal import resample, resample_poly
from fractions import Fraction

def ts_resample(timeseries, TR_new, TR_org=None, time=None, method='poly', axis=-1):

    n_tp = timeseries.shape[axis] # number of timepoints
    ratio = TR_org/TR_new

    if (TR_org == None) and (time == None):
        raise Exception('TR_org and time cannot be None at the same time')
    
    if TR_org == None:
        TR_org = time[1]-time[0]

    if method == 'fft':
        ts_resampled = resample(timeseries, t=time, axis=axis, num=int(TR_org/TR_new*n_tp))
    elif method == 'poly':
        fraq = Fraction(ratio)
        ts_resampled = resample_poly(timeseries, up=fraq.numerator, down=fraq.denominator, axis=axis)
    #elif method == '':

    time_resampled = np.arange(0, n_tp*TR_org, TR_new)
    return ts_resampled, time_resampled


# load surface file to resample
func_gii = nib.load(snakemake.input[0])
data = np.vstack([d.data for d in func_gii.darrays]).T

# get specified target and original TRs
TR_new = snakemake.params.target_TR
TR_org = snakemake.params.original_TR

# resample the functional data
data_resampled, _ = ts_resample(data, TR_new, TR_org)

# save to gifti (keeping the input structure, only modifying data)
# initiallize the new darrays
darrays = list()
for ts in range(data_resampled.shape[1]):
    # create a GIFTI data array and add to darrays
    darrays.append(nib.gifti.GiftiDataArray(data=np.squeeze(data_resampled[:,ts].astype('float32')))) # No intent code has specified (may want to use normal)
func_gii.darrays = darrays
# save GIFTI image to output
nib.save(func_gii, snakemake.output[0])

