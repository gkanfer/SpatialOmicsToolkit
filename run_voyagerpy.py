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


path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
FilePrefix = "_072824" 
andata = sc.read_h5ad(os.path.join(pathout, "andata_save.h5ad"))
andata.uns['config'] = OrderedDict()
andata.uns['config'] = OrderedDict()
andata.uns["config"]["secondary_var_names"] = andata.var_names
from scipy.sparse import csr_matrix
andata_sub = andata.copy()
andata_sub.X = csr_matrix(andata_sub.X)
andata_sub = sc.pp.subsample(andata_sub, n_obs=200_000,copy=True)
scale = 1
visium_spots = gpd.GeoSeries.from_xy(andata_sub.obsm['spatial'][:,0], andata_sub.obsm['spatial'][:,1]).scale(scale, scale, origin=(0, 0))
_ = vp.spatial.set_geometry(andata_sub, geom="spot_poly", values=visium_spots)

andata_sub.X = andata_sub.layers['log']
sc.pp.scale(andata_sub, max_value=10)
sc.pp.pca(andata_sub, n_comps=15)
sc.pp.neighbors(
    andata_sub,
    n_neighbors=25,
    n_pcs=15,
    use_rep='X_pca',
    knn=True,
    random_state=29403943,
    method='umap', # one of umap, gauss, rapids
    metric='cosine', # many metrics available,
    key_added='knn'
)
dist = andata_sub.obsp['knn_distances'].copy()
#dist.data[dist.data == 0] = 0.000001
dist.data = 1 / dist.data

# row normalize the matrix, this makes the matrix dense.
dist /= dist.sum(axis=1)

# convert dist back to sparse matrix
from scipy.sparse import csr_matrix
andata_sub.obsp["knn_weights"] = csr_matrix(dist)

del dist

knn_graph = "knn_weights"

# adata.obsp["knn_connectivities"] represent the edges, while adata.opsp["knn_weights"] represent the weights
andata_sub.obsp["knn_connectivities"] = (andata_sub.obsp[knn_graph] > 0).astype(int)
vp.spatial.set_default_graph(andata_sub, "knn_weights")
vp.spatial.to_spatial_weights(andata_sub, graph_name=knn_graph)
qc_features = ["total_counts"]
morans = vp.spatial.moran(andata_sub, qc_features, graph_name=knn_graph)
andata_sub.uns['spatial']['moran'][knn_graph].loc[qc_features, ["I"]]

qc_features = ["total_counts"]
vp.spatial.compute_spatial_lag(
    andata_sub,
    qc_features,
    graph_name=knn_graph,
    inplace=True
)

_ = vp.spatial.local_moran(andata_sub, qc_features, graph_name=knn_graph)
with PdfPages(os.path.join(pathout, 'barcode_hist_aoutocorlation.pdf')) as pdf:
    axs = vp.plt.plot_barcode_histogram(
    andata_sub,
    qc_features,
    obsm="local_moran",
    color_by='cluster',
    log=True,
    histtype='line',
    bins=10)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    pdf.savefig()
    plt.close()   