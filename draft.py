import cupy as cp
import cupyx
import scanpy as sc
import numpy as np
import pandas as pd
from cupyx.scipy.sparse import csr_matrix
import os
from PIL import Image
from sklearn.linear_model import LinearRegression
import pickle
import pickle
import esda
import pandas as pd
import geopandas as gpd
from geopandas import GeoDataFrame
import libpysal as lps
from libpysal.weights import W
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
import seaborn as sns 
from scipy.sparse import csr_matrix
import pickle

path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
FilePrefix = "_072824" 

grid = sc.read_h5ad(os.path.join(pathout, "grid_save.h5ad"))
file_path = os.path.join(pathout, "grid_uns_mtracies.pkl")

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
sparse_matrix = grid.X
row_sums = sparse_matrix.sum(axis=1)
grid.obs['n_counts'] = np.array(row_sums).flatten()

grid.layers['counts'] = grid.X.copy()
sc.pp.normalize_total(grid)
sc.pp.log1p(grid)
grid.layers['log'] = grid.X.copy()
sc.pp.scale(grid, max_value=10)
sc.tl.pca(grid, use_highly_variable=False, n_comps=15, random_state=1337)
sc.pp.neighbors(
    grid,
    n_neighbors=20,
    n_pcs=9,
    use_rep='X_pca',
    knn=True,
    random_state=29,
    method='gauss', # one of umap, gauss, rapids
    metric='euclidean', # many metrics available,
    key_added='knn'
)

dist = andata016_.obsp['knn_distances'].copy()
dist.data = 1 / dist.data

# row normalize the matrix, this makes the matrix dense.
dist /= dist.sum(axis=1)

# convert dist back to sparse matrix
from scipy.sparse import csr_matrix
andata016_.obsp["knn_weights"] = csr_matrix(dist)

del dist

knn_graph = "knn_weights"

# adata.obsp["knn_connectivities"] represent the edges, while adata.opsp["knn_weights"] represent the weights
andata016_.obsp["knn_connectivities"] = (andata016_.obsp[knn_graph] > 0).astype(int)
vp.spatial.set_default_graph(andata016_, knn_graph)
vp.spatial.to_spatial_weights(andata016_, graph_name=knn_graph)