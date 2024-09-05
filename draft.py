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
import libpysal as lps
import scipy.sparse as sp


########## Xenium ###################

# path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
# pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
# FilePrefix = "_072824" 
# andata = sc.read_h5ad(os.path.join(pathout, "andata_save.h5ad"))
# andata.uns['config'] = OrderedDict()
# andata.uns['config'] = OrderedDict()
# andata.uns["config"]["secondary_var_names"] = andata.var_names
# from scipy.sparse import csr_matrix
# andata_sub = andata.copy()
# andata_sub.X = csr_matrix(andata_sub.X)
# andata_sub = sc.pp.subsample(andata_sub, n_obs=7000,copy=True)
# scale = 1
# visium_spots = gpd.GeoSeries.from_xy(andata_sub.obsm['spatial'][:,0], andata_sub.obsm['spatial'][:,1]).scale(scale, scale, origin=(0, 0))
# _ = vp.spatial.set_geometry(andata_sub, geom="spot_poly", values=visium_spots)


# andata_sub.X = andata_sub.layers['log']
# sc.pp.scale(andata_sub, max_value=10)

# sc.pp.pca(andata_sub, n_comps=15)
# sc.pp.neighbors(
#     andata_sub,
#     n_neighbors=25,
#     n_pcs=15,
#     use_rep='X_pca',
#     knn=True,
#     random_state=29403943,
#     method='umap', # one of umap, gauss, rapids
#     metric='cosine', # many metrics available,
#     key_added='knn'
# )


# dist = andata.obsp['distances'].tocsr()
# #dist.data[dist.data == 0] = 0.000001
# dist.data = 1 / dist.data

# # row normalize the matrix, this makes the matrix dense.
# dist /= dist.sum(axis=1)

# # convert dist back to sparse matrix
# from scipy.sparse import csr_matrix
# andata_sub.obsp["knn_weights"] = csr_matrix(dist)

# del dist

# knn_graph = "knn_weights"

# # adata.obsp["knn_connectivities"] represent the edges, while adata.opsp["knn_weights"] represent the weights
# andata_sub.obsp["knn_connectivities"] = (andata_sub.obsp[knn_graph] > 0).astype(int)
# vp.spatial.set_default_graph(andata_sub, "knn_weights")
# vp.spatial.to_spatial_weights(andata_sub, graph_name=knn_graph)

################# Xinum Voyger coustom ###################

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
andata_sub = sc.pp.subsample(andata_sub, n_obs=7000,copy=True)
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
W = lps.weights.W.from_adjlist(graph_df_filtered)
W.set_transform("r")


# ########## Visium ###################

# import cupy as cp
# import cupyx
# import scanpy as sc
# import scanpy as sc
# import numpy as np
# from cupyx.scipy.sparse import csr_matrix
# import os
# from PIL import Image
# import matplotlib.pyplot as plt
# from matplotlib.colors import ListedColormap
# from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.backends.backend_pdf import PdfPages
# import seaborn as sns
# import random
# import pandas as pd
# import voyagerpy as vp
# import geopandas as gpd
# import libpysal as lps
# from collections import OrderedDict

# path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Human_Colon_Cancer_binned_outputs/binned_outputs/square_016um"
# pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
# andata = sc.read_h5ad(os.path.join(pathout, "andata_save_colon.h5ad"))
# andata.uns['config'] = OrderedDict()
# andata.uns['config'] = OrderedDict()
# andata.uns["config"]["secondary_var_names"] = andata.var_names
# def load_matrix(andata,pathout,npz_file = "obsp_distances.npz",mode = 'sparse', mat_name = 'distances'):
#     npzfile = cp.load(os.path.join(pathout, npz_file))
#     data = cp.array(npzfile['data'])
#     indices = cp.array(npzfile['indices'])
#     indptr = cp.array(npzfile['indptr'])
#     shape = tuple(npzfile['shape'])
    
#     # Reconstruct the sparse matrix
#     sparse_matrix_distances = csr_matrix((data, indices, indptr), shape=shape)
#     if mode== 'sparse':
#         andata.obsp[mat_name] = sparse_matrix_distances
#     else:
#         andata.obsp[mat_name] = sparse_matrix_distances.get()
#     return andata
# andata = load_matrix(andata,pathout,npz_file = "obsp_distances_large_colon.npz",mode = 'sparse', mat_name = 'distances')
# andata = load_matrix(andata,pathout,npz_file = "obsp_connectivities_large_colon.npz",mode = 'sparse', mat_name = 'connectivities')

# import scipy.sparse as sp
# sparse_dist_matrix = andata.obsp['distances'].tocsr()
# sparse_inv_matrix = sparse_dist_matrix.copy()
# sparse_inv_matrix.data = 1 / sparse_inv_matrix.data
# sparse_inv_matrix.data[sparse_inv_matrix.data == float('inf')] = 0

# import numpy as np
# import pandas as pd

# # Convert the sparse matrix to COOrdinate format
# sparse_inv_matrix_coo = sparse_inv_matrix.tocoo()

# # Extract the row (focal) and column (neighbor) indices of non-zero entries
# focal = sparse_inv_matrix_coo.row
# neighbors = sparse_inv_matrix_coo.col

# focal = focal.get()
# neighbors = neighbors.get() 

# idx = np.array(andata.obs_names,dtype=str) # Assuming this is a pandas Index or a list-like structure

# # Convert sparse matrix values to a 1D array explicitly
# weights = sparse_inv_matrix_coo.data
# weights = weights.get()



# # Create a DataFrame with focal, neighbor, and weight information
# graph_df = pd.DataFrame({
#     "focal": idx[focal],
#     "neighbor": idx[neighbors],
#     "weight": weights  # The actual non-zero values (inverted distances)
# })

# # Display the DataFrame to check
# graph_df.head()

# graph_df_filtered = graph_df[graph_df['weight'] != 0]

# W = lps.weights.W.from_adjlist(graph_df_filtered)
# W.set_transform("r")
