import debugpy
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import os
import gzip
import numpy as np
import cell2location

# The model and output h5ad can be loaded later like this:
results_folder = './out/cell2loc/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

def calcQCmat(andata):
    andata.var_names_make_unique()
    andata.var["mt"] = andata.var_names.str.startswith("mt-")
    andata.var["ribo"] = andata.var_names.str.startswith(("RPS", "RPL"))
    andata.var["hb"] = andata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(andata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    return andata
path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
andata016_ = sc.read_visium(path=path_016)
andata016 = calcQCmat(andata016_)
andata016 = andata016[:,andata016.var.n_cells_by_counts > 50]
andata016.raw = andata016.copy()
sc.pp.normalize_total(andata016)
sc.pp.log1p(andata016)
log1p_data = andata016.X.todense()
sc.pp.highly_variable_genes(andata016)
sc.pp.scale(andata016)
andata016.obsm['spatial'] = np.array(andata016.obsm['spatial'], dtype=np.float64)
sc.pp.pca(andata016, n_comps=10)
sc.pp.neighbors(andata016)
sc.tl.leiden(andata016, key_added="clusters" ,resolution=0.7 , directed=False, n_iterations=2)


andata016.X = andata016.raw.X
andata016.var.set_index('gene_ids', drop=True, inplace=True)

intersect = np.intersect1d(andata016.var_names, inf_aver.index)
adata_vis = andata016[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
import pickle
path_out = './pymc_model'
df = pd.DataFrame(data=adata_vis.X.todense(), columns=adata_vis.var.index.tolist())
with open(os.path.join(path_out, 'brain_16mm_vshd.pkl'), 'wb') as buff:
            pickle.dump({'dsg':df,'wsf':inf_aver}, buff)