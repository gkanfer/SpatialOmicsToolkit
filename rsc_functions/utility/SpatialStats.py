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
import cupyx

def compute_spatial_lag(andata,feature,knn_key = 'knn_'):
    '''
        rsc.pp.neighbors(andata, n_pcs=15, use_rep='X_pca', n_neighbors=30,key_added = 'knn')
    '''
    dist = andata.obsp[f'{knn_key}distances'].copy()
    epsilon = 1e-10
    dist.data = 1 / (dist.data + epsilon)
    # row normalize the matrix, this makes the matrix dense.
    dist /= dist.sum(axis=1)
    
    andata.obsp["knn_weights"] = cupyx.scipy.sparse.csr_matrix(dist)
    del dist
    knn_graph = "knn_weights"
    andata.obsp[f"{knn_key}connectivities"] = (andata.obsp[knn_graph] > 0).astype(cp.float32)    
    dists = andata.obsp[f'{knn_key}connectivities']
    dists = cupyx.scipy.sparse.csr_matrix.toarray(dists)
    lagged_feat = f"lagged_{feature}"
    x = cp.array(andata.obs[feature].values)
    # Perform the matrix multiplication in CuPy
    lagged_total_counts = dists.dot(x)
    lagged_total_counts_numpy = lagged_total_counts.get()
    andata.obs[lagged_feat] = lagged_total_counts_numpy
    del dists
    return andata
    

