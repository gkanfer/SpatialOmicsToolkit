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
import libpysal as lps
from collections import OrderedDict
import scipy.sparse as sp
import pickle

path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Human_Colon_Cancer_binned_outputs/binned_outputs/square_016um"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
andata = sc.read_h5ad(os.path.join(pathout, "andata_save_colon.h5ad"))
andata.uns['config'] = OrderedDict()
andata.uns['config'] = OrderedDict()
andata.uns["config"]["secondary_var_names"] = andata.var_names
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
andata = load_matrix(andata,pathout,npz_file = "obsp_distances_large_colon.npz",mode = 'cupyx', mat_name = 'distances')
andata = load_matrix(andata,pathout,npz_file = "obsp_connectivities_large_colon.npz",mode = 'cupyx', mat_name = 'connectivities')
andata

with open(os.path.join(pathout,"andata_uns_mtracies__colon.pkl"), 'rb') as buff:
    andata_uns_mtracies__colon = pickle.load(buff)
andata.uns['rank_genes_groups'] = andata_uns_mtracies__colon['rank_genes_groups']
del andata_uns_mtracies__colon
with open(os.path.join(pathout,"stlearn_uns_mtracies_colon.pkl"), 'rb') as buff:
    stlearn_uns_mtracies_colon = pickle.load(buff)
andata.uns['lrfeatures'] = stlearn_uns_mtracies_colon['lrfeatures']
andata.uns['lr_summary'] = stlearn_uns_mtracies_colon['lr_summary']
del stlearn_uns_mtracies_colon
with open(os.path.join(pathout,"stlearn_obsm_mtracies_colon.pkl"), 'rb') as buff:
    stlearn_obsm_mtracies_colon = pickle.load(buff)
andata.obsm['spot_neighbours'] = stlearn_obsm_mtracies_colon['spot_neighbours']
andata.obsm['spot_neigh_bcs'] = stlearn_obsm_mtracies_colon['spot_neigh_bcs']
andata.obsm['lr_scores'] = stlearn_obsm_mtracies_colon['lr_scores']
andata.obsm['p_vals'] = stlearn_obsm_mtracies_colon['spot_neigh_bcs']
andata.obsm['-log10(p_adjs)'] = stlearn_obsm_mtracies_colon['-log10(p_adjs)']
andata.obsm['lr_sig_scores'] = stlearn_obsm_mtracies_colon['lr_sig_scores']
del stlearn_obsm_mtracies_colon

from scipy.sparse import csr_matrix
andata_sub = andata.copy()
andata_sub.X = csr_matrix(andata_sub.X)
andata_sub = sc.pp.subsample(andata_sub, n_obs=100_000,copy=True)


from cupyx.scipy.sparse import csr_matrix
andata_sub.obsp['distances'] = csr_matrix(andata_sub.obsp['distances'])
andata_sub.obsp['distances'] = csr_matrix(andata_sub.obsp['connectivities'])

import scipy.sparse as sp
sparse_dist_matrix = andata_sub.obsp['distances'].tocsr()
sparse_inv_matrix = sparse_dist_matrix.copy()
sparse_inv_matrix.data = 1 / sparse_inv_matrix.data
sparse_inv_matrix.data[sparse_inv_matrix.data == float('inf')] = 0

# Convert the sparse matrix to COOrdinate format
sparse_inv_matrix_coo = sparse_inv_matrix.tocoo()

# Extract the row (focal) and column (neighbor) indices of non-zero entries
focal = sparse_inv_matrix_coo.row
neighbors = sparse_inv_matrix_coo.col

focal = focal.get()
neighbors = neighbors.get() 

idx = np.array(andata_sub.obs_names,dtype=str) # Assuming this is a pandas Index or a list-like structure

# Convert sparse matrix values to a 1D array explicitly
weights = sparse_inv_matrix_coo.data
weights = weights.get()


# Create a DataFrame with focal, neighbor, and weight information
graph_df = pd.DataFrame({
    "focal": idx[focal],
    "neighbor": idx[neighbors],
    "weight": weights  # The actual non-zero values (inverted distances)
})

# Display the DataFrame to check
graph_df.head()

graph_df_filtered = graph_df[graph_df['weight'] != 0]

neighbor_counts = graph_df_filtered.groupby("focal")["neighbor"].count()
# Identify connected nodes (nodes with at least 2 neighbor)
connected_nodes = neighbor_counts[neighbor_counts > 3].index
# Filter graph_df to keep only connected nodes
graph_df_filtered = graph_df_filtered[graph_df_filtered['focal'].isin(connected_nodes)]

W = lps.weights.W.from_adjlist(graph_df_filtered)
W.set_transform("r")

knn_graph = "knn_weights"
andata_sub.obsp["knn_weights"] = sparse_inv_matrix.get()

andata_sub.uns.setdefault("spatial", {})
andata_sub.uns["spatial"][knn_graph] = W

qc_features = ["total_counts"]
morans = vp.spatial.moran(andata_sub, qc_features, graph_name=knn_graph)
andata_sub.uns['spatial']['moran'][knn_graph].loc[qc_features, ["I"]]

ylag = lps.weights.lag_spatial(W, andata_sub.obs['total_counts'].values)
andata_sub.obs['lagged_total_counts'] = lps.weights.lag_spatial(W, andata_sub.obs['total_counts'].values)
_ = vp.spatial.local_moran(andata_sub, qc_features, graph_name=knn_graph)

hvg = andata_sub.var[andata_sub.var['highly_variable'].values].index
vp.spatial.moran(andata_sub, feature=hvg, dim='var', graph_name=knn_graph)

hvgs_moransI = andata_sub.uns['spatial']['moran'][knn_graph].loc[hvg, 'I']
andata_sub.var.loc[hvg, "moran"] = hvgs_moransI
andata_sub.var.loc[:, ["moran"]] = np.nan_to_num(andata_sub.var.loc[:, ["moran"]],0.0)
mat_names = np.ravel(pd.DataFrame(andata_sub.uns['rank_genes_groups']['names']).values)
marker_genes = mat_names[:20].tolist()
andata_sub.var['symbol'] = andata_sub.var['gene_ids'].values

with PdfPages(os.path.join(pathout, 'plot_features_histogram_marker_genes_colon.pdf')) as pdf:
    _ = vp.plt.plot_features_histogram(
        andata_sub,
        "moran",
        bins=50,
        log=False,
        histtype="bar",
        markers=marker_genes,
        show_symbol=False)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    pdf.savefig()
    plt.close()
    
cluster_num = len(np.unique(andata_sub.obs['cluster'].values))
_ = vp.spatial.local_moran(andata_sub, marker_genes, graph_name=knn_graph)
marker_genes_symbols = andata_sub.var.loc[marker_genes, "symbol"].tolist()

with PdfPages(os.path.join(pathout, 'plot_barcode_histogram_marker_genes_colon.pdf')) as pdf:
    _ = vp.plt.plot_barcode_histogram(
        andata_sub,
        marker_genes,
        color_by='cluster',
        obsm='local_moran',
        histtype='line',
        figsize=(15,15),
        subplot_kwargs=dict(layout='constrained'),
        label=marker_genes,
        ncol=2
    )
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    pdf.savefig()
    plt.close()
    