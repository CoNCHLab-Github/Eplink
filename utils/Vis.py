from surfplot import Plot
import matplotlib.pyplot as plt
import os

def plot_maps(map, atlas=None, mesh='very_inflated', surf_mesh=None, range=None, layout='grid' , cmap='viridis', cbar=True, title=None, clabel=None, ax=None):
    # Surface mesh files
    if surf_mesh is None:
        surf_mesh = os.path.join('ISC-pipeline','resources','fsLR_32K_surfaces',f'fs_LR.32k.{{hemi}}.{mesh}.surf.gii')

    cbar_kws = {'location': 'bottom', 'label_direction': 0, 'decimals': 2,
    'fontsize': 8, 'n_ticks': 7, 'shrink': 0.5, 'aspect': 40}
    
    # plot surface mesh
    p = Plot(surf_lh=surf_mesh.format(hemi='L'), surf_rh=surf_mesh.format(hemi='R'), brightness = 0.7, size=(600, 600), zoom=1.2, layout=layout)
    # plot the map
    p.add_layer({'left': map.format(hemi='L'), 'right': map.format(hemi='R')}, cbar=cbar, cbar_label = clabel, color_range=range, cmap=cmap)
    # plot atlas borders
    if atlas != None:
        p.add_layer({'left': atlas.format(hemi='L'), 'right': atlas.format(hemi='R')}, cmap='gray', as_outline=True, cbar=False)

    if ax is None:
        fig = p.build(cbar_kws=cbar_kws)
        fig.suptitle(title);
    else:
        pr = p.render()
        pr._check_offscreen()
        x = pr.to_numpy(transparent_bg=True, scale=(2,2))

        ax.imshow(x)
        ax.axis('off')
        ax.set_title(title)
        if cbar:
            # cbar_kws = {} if cbar_kws is None else cbar_kws
            plt.sca(ax)
            p._add_colorbars(**cbar_kws)
    return p
