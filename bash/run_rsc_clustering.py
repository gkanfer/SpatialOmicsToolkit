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
import cupy as cp
import os
import time
import rapids_singlecell as rsc
import numpy as np
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
from rsc_functions.reports.plot import plot_spatial

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
import os

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
rsc.pp.scale(andata)

rsc.pp.highly_variable_genes(andata)
rsc.pp.pca(andata, n_comps=15,random_state=1337,mask_var="highly_variable")
rsc.pp.neighbors(andata,n_pcs=15,use_rep='X_pca',n_neighbors=15)
rsc.tl.leiden(andata,random_state=1337,resolution=1,key_added='cluster') 

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(os.path.join(pathout, f'Report_spatialPlot_test.pdf')) as pdf:
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    plot_spatial(andata,ax = ax, features = None, title = '', xlab = '',ylab ='',size = 2)
    fig.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    pdf.savefig()
    plt.close()

