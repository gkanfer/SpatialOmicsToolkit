import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import os
import gzip
import numpy as np
import scanpy as sc
import squidpy as sq
import cupy as cp
import cupyx
import os
import time
import rapids_singlecell as rsc
import numpy as np
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
import cupy

rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)
import zarr
from collections import OrderedDict
from scipy.sparse import csr_matrix
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from scipy.sparse import csr_matrix
import scipy
import anndata
from collections import OrderedDict
from rsc_functions.utility.applyqc import applyqc
from rsc_functions.reports.plot import plot_spatial,plot_spatial_data
from rsc_functions.utility.rank_genes_groups import return_markers,rank_genes_groups
from rsc_functions.reports.plot import plot_expression


path_xenium = os.path.join(path,"cell_feature_matrix.h5")
path_cells = os.path.join(path,"cells.zarr.zip")
adata = sc.read_10x_h5(path_xenium)
rsc.get.anndata_to_GPU(adata)
rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="mt-")
rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT"])
def open_zarr(path: str) -> zarr.Group:
    store = (zarr.ZipStore(path, mode="r") if path.endswith(".zip") else zarr.DirectoryStore(path))
    return zarr.group(store=store)
root = open_zarr(path_cells)
column_names = dict(root['cell_summary'].attrs.items())['column_names']
def build_obs(andata,root,column_names):
    for i in range(len(column_names)):
        andata.obs[str(column_names[i])] = np.array(root["cell_summary"])[:,i]
    spatial = andata.obs[["cell_centroid_x", "cell_centroid_y"]]
    adata.obsm["spatial"] = spatial.values
    return andata
andata = build_obs(adata,root,column_names)
andata.var_names_make_unique()
andata.obsm['spatial'] = np.array(andata.obsm['spatial'], dtype=np.float64)
andata.uns['config'] = OrderedDict()
andata.uns["config"]["secondary_var_names"] = andata.var_names
rsc.pp.flag_gene_family(andata, gene_family_name="MT", gene_family_prefix="mt-")
rsc.pp.calculate_qc_metrics(andata, qc_vars=["MT"])
rsc.pp.filter_cells(andata, min_count=10,qc_var = 'total_counts')
rsc.pp.filter_genes(andata, min_count=5)
andata.layers['counts'] = andata.X.copy()
rsc.pp.normalize_total(andata)
rsc.pp.log1p(andata)
andata.layers['log'] = andata.X.copy()
rsc.pp.highly_variable_genes(andata, n_top_genes=1500, flavor="seurat_v3", layer="log")
andata = andata[:, andata.var["highly_variable"]]
rsc.pp.scale(andata, max_value=10)
rsc.pp.pca(andata, n_comps=30,random_state=1337, use_highly_variable=False)


rsc.pp.neighbors(andata, n_pcs=15, use_rep='X_pca', n_neighbors=30)
rsc.tl.leiden(andata, random_state=1337, resolution=0.5, key_added='cluster') 


####################################################################################################
################################## write to andata #################################################
andata_save = andata.copy()
andata_save.X = andata_save.layers['log']
del andata_save.uns
del andata_save.obsm
del andata_save.varm
del andata_save.layers
del andata_save.obsp
andata_save.write_h5ad(os.path.join(pathout, "andata_save.h5ad"))

####################################################################################################
################################## read #################################################

andata = sc.read_h5ad(os.path.join(pathout, "andata_save.h5ad"))



####################################################################################################
################################## write knn matrix using gpu ######################################

import cupy as cp
import os

# Assuming `andata.obsp['distances']` is a CuPy sparse matrix in CSR format
sparse_matrix = andata.obsp['distances']

# Save the sparse matrix components
cp.savez(os.path.join(pathout, "obsp_distances.npz"), 
         data=sparse_matrix.data, 
         indices=sparse_matrix.indices, 
         indptr=sparse_matrix.indptr, 
         shape=sparse_matrix.shape)

####################################################################################################
################################## read knn matrix using gpu ######################################


import cupy as cp
from cupyx.scipy.sparse import csr_matrix

# Load the saved sparse matrix components
npzfile = cp.load(os.path.join(pathout, "obsp_distances.npz"))
data = cp.array(npzfile['data'])
indices = cp.array(npzfile['indices'])
indptr = cp.array(npzfile['indptr'])
shape = tuple(npzfile['shape'])

# Reconstruct the sparse matrix
sparse_matrix_gpu = csr_matrix((data, indices, indptr), shape=shape)

# Verify the result
print(sparse_matrix_gpu)
print(f"Shape: {sparse_matrix_gpu.shape}")

'''
(0, 0)	0.0
  (0, 383156)	2.3886120319366455
  (0, 356986)	2.5408105850219727
  (0, 56709)	2.7588164806365967
  (0, 681883)	2.803593635559082
  (0, 372891)	2.880432367324829
  (0, 56023)	2.941385507583618
  (0, 281931)	2.9922263622283936
  (0, 441360)	3.016859292984009
  (0, 618469)	3.1007328033447266
  (0, 138264)	3.1349611282348633
  (0, 300878)	3.144524097442627
  (0, 604064)	3.204695463180542
  (0, 388517)	3.2333505153656006
  
  (705297, 511425)	2.4162704944610596
  (705297, 473073)	2.5159220695495605
  (705297, 420205)	2.5167438983917236
  (705297, 407297)	2.5952112674713135
  (705297, 447754)	2.6227433681488037
  (705297, 475440)	2.6888139247894287
  (705297, 464941)	2.702136993408203
  (705297, 635485)	2.706777572631836
  (705297, 466814)	2.746880054473877
  (705297, 519314)	2.795445203781128
  (705297, 702324)	2.803297281265259
  (705297, 702608)	2.8066372871398926
  (705297, 404978)	2.8075151443481445
  (705297, 662373)	2.8405449390411377
  (705297, 632587)	2.841468572616577
  (705297, 566602)	2.8506875038146973
  (705297, 622078)	2.860783338546753
  (705297, 456247)	2.8686158657073975
  (705297, 637393)	2.876199245452881
  (705297, 580266)	2.907470464706421
  (705297, 658675)	2.910426378250122
  (705297, 580265)	2.925283193588257
  (705297, 512490)	2.9255359172821045
  (705297, 663562)	2.938891649246216
  (705297, 410359)	2.9491074085235596
Shape: (705298, 705298)
'''
