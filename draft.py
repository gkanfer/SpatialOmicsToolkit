import cupy as cp
import cupyx
import scanpy as sc
import scanpy as sc
import numpy as np
from cupyx.scipy.sparse import csr_matrix
import os
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import random
import pandas as pd
import voyagerpy as vp
import geopandas as gpd
from collections import OrderedDict
import pickle

path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Human_Colon_Cancer_binned_outputs/binned_outputs/square_016um"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
grid = sc.read_h5ad(os.path.join(pathout, "grid_save_colon.h5ad"))
file_path = os.path.join(pathout, "grid_uns_mtracies_colon.pkl")
if os.path.getsize(file_path) > 0:
    with open(file_path, 'rb') as buff:
        grid_uns_mtracies = pickle.load(buff)
else:
    print("File is empty. Cannot load data.")
    grid_uns_mtracies = None
grid.uns = {}
grid.uns['cluster'] = pd.read_csv(os.path.join(pathout, "cluster.csv"))
grid.uns['grid_counts'] = grid_uns_mtracies['grid_counts']
grid.uns['grid_xedges'] = grid_uns_mtracies['grid_xedges']
grid.uns['grid_yedges'] = grid_uns_mtracies['grid_yedges']
grid.uns['lrfeatures'] = pd.read_csv(os.path.join(pathout, "lrfeatures.csv"))
grid.uns['lr_summary'] = pd.read_csv(os.path.join(pathout, "lr_summary.csv"))
#grid.X = grid.layers['count']
sparse_matrix = grid.X
row_sums = sparse_matrix.sum(axis=1)
grid.obs['n_counts'] = np.array(row_sums).flatten()
grid.layers['counts'] = grid.X.copy()
sc.pp.normalize_total(grid)
sc.pp.log1p(grid)
grid.layers['log'] = grid.X.copy()
sc.pp.scale(grid, max_value=10)

scale = 1
visium_spots = gpd.GeoSeries.from_xy(grid.obsm['spatial'][:,0], grid.obsm['spatial'][:,1]).scale(scale, scale, origin=(0, 0))
_ = vp.spatial.set_geometry(grid, geom="spot_poly", values=visium_spots)
grid.uns['config'] = OrderedDict()
grid.uns["config"]["secondary_var_names"] = grid.var_names
sc.pp.pca(grid, n_comps=15)
sc.pp.neighbors(
    grid,
    n_neighbors=25,
    n_pcs=15,
    use_rep='X_pca',
    knn=True,
    random_state=29403943,
    method='umap', # one of umap, gauss, rapids
    metric='cosine', # many metrics available,
    key_added='knn'
)
dist = grid.obsp['knn_distances'].copy()
dist.data = 1 / dist.data

# row normalize the matrix, this makes the matrix dense.
dist /= dist.sum(axis=1)

# convert dist back to sparse matrix
from scipy.sparse import csr_matrix
grid.obsp["knn_weights"] = csr_matrix(dist)

del dist

knn_graph = "knn_weights"

# adata.obsp["knn_connectivities"] represent the edges, while adata.opsp["knn_weights"] represent the weights
grid.obsp["knn_connectivities"] = (grid.obsp[knn_graph] > 0).astype(int)
vp.spatial.set_default_graph(grid, knn_graph)
vp.spatial.to_spatial_weights(grid, graph_name=knn_graph)