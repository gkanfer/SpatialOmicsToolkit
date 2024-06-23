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
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

# The model and output h5ad can be loaded later like this:
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_ref = sc.read_h5ad(adata_file)
# mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)