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
path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
FilePrefix = "_072824" 
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

sc.pp.subsample(andata,copy=False,n_obs=100_000)

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

#rsc.tl.rank_genes_groups_logreg(andata, groupby="cluster",groups='all')

# andata = sc.read_h5ad(os.path.join(pathout, "andata_save.h5ad"))
# # data is looged
rsc.tl.rank_genes_groups_logreg(andata, groupby="cluster", use_raw=False)
# sc.tl.rank_genes_groups(andata, groupby="cluster", method="logreg")
rsc.get.anndata_to_CPU(andata)
sc.pl.rank_genes_groups_dotplot(andata, groupby="cluster",n_genes = 3)