from skimage import io
import numpy as np
import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import gzip
import numpy as np
import celltypist
from celltypist import models


plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

def calcQCmat(andata):
    andata.var_names_make_unique()
    andata.var["mt"] = andata.var_names.str.startswith("mt-")
    andata.var["ribo"] = andata.var_names.str.startswith(("RPS", "RPL"))
    andata.var["hb"] = andata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(andata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    return andata
path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out"
andata016_ = sc.read_visium(path=path_016)
andata016 = calcQCmat(andata016_)
print(f"{np.shape(andata016.X.todense())}")
sc.pp.filter_cells(andata016, min_counts = 50)
sc.pp.filter_cells(andata016, min_genes = 80)

andata016 = andata016[:,andata016.var.n_cells_by_counts > 50]
andata016 = andata016[andata016.obs["pct_counts_mt"] < 20]

sc.pp.normalize_total(andata016)
sc.pp.log1p(andata016)
log1p_data = andata016.X.todense()
sc.pp.highly_variable_genes(andata016)
sc.pp.scale(andata016)
andata016.obsm['spatial'] = np.array(andata016.obsm['spatial'], dtype=np.float64)
sc.pp.pca(andata016, n_comps=10)
sc.pp.neighbors(andata016)
sc.tl.umap(andata016)
sc.tl.leiden(andata016, key_added=f'clusters', flavor="igraph", directed=False, resolution=0.5, n_iterations=2)

sc.tl.rank_genes_groups(
    andata016, groupby="clusters", method="wilcoxon", key_added="dea_clusters"
)
