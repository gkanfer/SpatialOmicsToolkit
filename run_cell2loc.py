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
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

outPath =  "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"

results_folder = './out/cell2loc/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

adata_ref = sc.read('38795736-71ed-4b5d-b7a1-91a9b60d52a2.h5ad')
adata_ref.X = adata_ref.raw.X

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()

cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Batch_ID',
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['assay']
                       )
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, use_gpu=True)

with PdfPages(os.path.join(outPath, f'ELBO_loss_history.pdf')) as pdf:
    mod.plot_history(20)
    
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000,'batch_size': 2500, 'use_gpu': True}
)

with PdfPages(os.path.join(outPath, f'adata_ref_QC_plots.pdf')) as pdf:
    mod.plot_QC()

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_ref.write(adata_file)
# adata_file

# The model and output h5ad can be loaded later like this:
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_ref = sc.read_h5ad(adata_file)
# mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
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


cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=20,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=200
)


mod.train(max_epochs=3000,
          # train using full data (batch_size=None)
          batch_size=3000,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

with PdfPages(os.path.join(outPath, f'cell2loc_model.pdf')) as pdf:
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])

# sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf', n_neighbors = 15)
# sc.tl.leiden(adata_vis, resolution=0.5)
# adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")
adata_vis.write(os.path.join(outPath, "adata_vis.h5ad"))



