import cupy as cp
import cupyx
import scanpy as sc
import stlearn as st
import scanpy as sc
import numpy as np
from cupyx.scipy.sparse import csr_matrix
import os
from PIL import Image

path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
FilePrefix = "_072824" 
andata = sc.read_h5ad(os.path.join(pathout, "andata_save.h5ad"))

def load_matrix(andata,pathout,npz_file = "obsp_distances.npz",mode = 'sparse', mat_name = 'distances'):
    npzfile = cp.load(os.path.join(pathout, npz_file))
    data = cp.array(npzfile['data'])
    indices = cp.array(npzfile['indices'])
    indptr = cp.array(npzfile['indptr'])
    shape = tuple(npzfile['shape'])
    
    # Reconstruct the sparse matrix
    sparse_matrix_distances = csr_matrix((data, indices, indptr), shape=shape)
    if mode== 'sparse':
        andata.obsp[mat_name] = sparse_matrix_distances
    else:
        andata.obsp[mat_name] = sparse_matrix_distances.get()
    return andata
andata = load_matrix(andata,pathout,npz_file = "obsp_connectivities_large.npz",mode = 'sparse', mat_name = 'distances')
andata = load_matrix(andata,pathout,npz_file = "obsp_distances_large.npz",mode = 'sparse', mat_name = 'connectivities')
### Calculating the number of grid spots we will generate

andata.obs["imagecol"] = andata.obs["cell_centroid_x"]
andata.obs["imagerow"] = andata.obs["cell_centroid_y"]
andata.uns["spatial"] = andata.obsm["spatial"]
def make_uns_spatial(adata):
    max_size = np.max([adata.obs["imagecol"].max(), adata.obs["imagerow"].max()])
    max_size = int(max_size + 0.1 * max_size)
    image = Image.new("RGBA", (max_size, max_size), (255, 255, 255, 255))
    imgarr = np.array(image)
    adata.uns["spatial"] = {}
    adata.uns["spatial"]["id1"] = {}
    adata.uns["spatial"]["id1"]["images"] = {}
    adata.uns["spatial"]["id1"]["images"]["lowres"] = imgarr
    adata.uns["spatial"]["id1"]["use_quality"] = "lowres"
    adata.uns["spatial"]["id1"]["scalefactors"] = {}
    adata.uns["spatial"]["id1"]["scalefactors"]["tissue_low_scalef" ] = 1
    adata.uns["spatial"]["id1"]["scalefactors"]["spot_diameter_fullres"] = 15
    adata.uns["spatial"]["id1"]["scalefactors"][
            "tissue_" + "lowres" + "_scalef"
        ] = 1
    return adata
andata = make_uns_spatial(andata)

n_ = 250
print(f'{n_} by {n_} has this many spots:\n', n_*n_)

### Gridding.
grid = st.tl.cci.grid(andata, n_row=n_, n_col=n_, use_label = 'cluster')
print( grid.shape ) # Slightly less than the above calculation, since we filter out spots with 0 cells.

# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(grid, lrs,
                min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                n_pairs=10_000, # Number of random pairs to generate; low as example, recommend ~10,000
                n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
                )

grid_save = grid.copy()
grid_save.uns['lr_summary'].to_csv(os.path.join(pathout, "lr_summary.csv"))
grid_save.uns['cluster'].to_csv(os.path.join(pathout, "cluster.csv"))
grid_save.uns['lrfeatures'].to_csv(os.path.join(pathout, "lrfeatures.csv"))

import pickle
with open(os.path.join(pathout,"grid_uns_mtracies.pkl"), 'wb') as buff:
     pickle.dump({'grid_counts':grid.uns['grid_counts'], 'grid_xedges':grid.uns['grid_xedges'], 'grid_yedges':grid.uns['grid_yedges']}, buff)

del grid_save.uns
grid_save.write_h5ad(os.path.join(pathout, "grid_save.h5ad"))

