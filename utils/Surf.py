import nibabel as nib
import numpy as np
import os

def roi2gii(roi_data, atlasfile, fpath, fname, mask=None):
    """Create GIFTI Image from ROI values and save it in the results."""
    # If data is not formatted as a dictionary we assume the first half corresponds to the left hemisphere
    if ~isinstance(roi_data,dict):  
        n_units = roi_data.shape[0]
        roi_data = {'L': roi_data[:n_units//2], 'R': roi_data[n_units//2:]}

    # Iterate on Left and Right hemispheres
    surface_values = dict()
    for hemi in ['L', 'R']:
        # Load atlas file for the hemisphere
        label_gii = nib.load(atlasfile.format(hemi=hemi))
        label_data = label_gii.darrays[0].data

        # idx2label = np.sort(np.unique(label_data))
        idx2label = list(label_gii.labeltable.get_labels_as_dict().keys())
        # Copy ROI value to vertices 
        n_vertices = len(label_data)
        surface_values[hemi] = np.zeros(n_vertices, dtype='float32')
        for roi_idx, value in enumerate(roi_data[hemi]):
            if (mask != None) and (not mask[hemi][roi_idx]):
                continue
            vertices = (label_data == idx2label[roi_idx])
            surface_values[hemi][vertices] = value
        
    vertex2gii(surface_values, fpath, fname)

    return surface_values

def vertex2gii(vertex_data, fpath, fname):
    """Create GIFTI Image from ROI values and save it in the results."""
    if not os.path.isdir(fpath):
        os.mkdir(fpath)
    # Iterate on Left and Right hemispheres
    for hemi in ['L', 'R']:
        # Create a GIFTI data array for the map
        gii_darray = nib.gifti.GiftiDataArray(data=vertex_data[hemi].astype('float32'), intent=nib.nifti1.intent_codes['NIFTI_INTENT_SHAPE'])
        # Create a GIFTI image for the map
        gii_image = nib.gifti.GiftiImage(darrays=[gii_darray])
        # Save GIFTI image
        nib.save(gii_image, os.path.join(fpath,f'{fname}_{hemi}.shape.gii'))

    