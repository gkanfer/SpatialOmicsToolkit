import scanpy as sc
import cupy as cp
import time
import rapids_singlecell as rsc
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)
import zarr
from collections import OrderedDict
from scipy.sparse import csr_matrix
from leidenalg import ModularityVertexPartition
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

adata = sc.read_10x_h5("/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp/cell_feature_matrix.h5")
adata.var_names_make_unique()
rsc.get.anndata_to_GPU(adata)
rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="mt-")
rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT"])

# Function to open a Zarr file
def open_zarr(path: str) -> zarr.Group:
    store = (zarr.ZipStore(path, mode="r") if path.endswith(".zip") else zarr.DirectoryStore(path))
    return zarr.group(store=store)
path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp/cells.zarr.zip"
root = open_zarr(path)
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

#rsc.pp.filter_cells(andata, min_count=10)
rsc.pp.filter_genes(andata, min_count=10)
andata.uns['config'] = OrderedDict()
andata.uns["config"]["secondary_var_names"] = andata.var_names
andata.layers['counts'] = andata.X.copy()
rsc.pp.normalize_total(andata)
rsc.pp.log1p(andata)
rsc.pp.highly_variable_genes(andata)
rsc.pp.pca(andata, n_comps=15,random_state=1337)
rsc.pp.neighbors(andata,n_pcs=15,use_rep='X_pca',n_neighbors=45)
rsc.tl.leiden(andata,random_state=1337,resolution=1,key_added='cluster')
rsc.tl.rank_genes_groups_logreg(andata, groupby="cluster",groups='all')
def z_to_p(z):
    return 2 * (1 - stats.norm.cdf(abs(z)))
p_values = pd.DataFrame.from_records(andata.uns['rank_genes_groups']['scores']).applymap(z_to_p)
p_values_flat = p_values.values.flatten()
_, pvals_corrected, _, _ = multipletests(p_values_flat, alpha=0.05, method='fdr_bh')
pvals_corrected_df = pd.DataFrame(pvals_corrected.reshape(p_values.shape), columns=p_values.columns, index=p_values.index)