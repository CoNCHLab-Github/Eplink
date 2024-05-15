from surfplot import Plot
import os

def plot_maps(map, atlas=None, mesh='very_inflated', range=None, layout='grid' , cmap='viridis', cbar=True, title=None):
    # Surface mesh files
    surf_mesh = os.path.join('ISC-pipeline','resources','fsLR_32K_surfaces',f'fs_LR.32k.{{hemi}}.{mesh}.surf.gii')

    kws = {'location': 'bottom', 'label_direction': 45, 'decimals': 2,
    'fontsize': 8, 'n_ticks': 7, 'shrink': 0.5, 'aspect': 40}
    
    # plot surface mesh
    p = Plot(surf_lh=surf_mesh.format(hemi='L'), surf_rh=surf_mesh.format(hemi='R'), brightness = 0.7, size=(1600, 400), zoom=1.2, layout=layout)
    # plot the map
    p.add_layer({'left': map.format(hemi='L'), 'right': map.format(hemi='R')}, cbar=cbar, color_range=range, cmap=cmap)
    # plot atlas borders
    if atlas != None:
        p.add_layer({'left': atlas.format(hemi='L'), 'right': atlas.format(hemi='R')}, cmap='gray', as_outline=True, cbar=False)

    fig = p.build(cbar_kws=kws)
    fig.suptitle(title);