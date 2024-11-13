import os
import tempfile

import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
# from scib_metrics.benchmark import Benchmarker

scvi.settings.seed = 0
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_4"
andata = sc.read_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_spatialleiden_harmony_scVI_.h5ad"))

scvi.model.SCVI.setup_anndata(andata, layer="counts", batch_key="batch")
model_scvi = scvi.model.SCVI(andata)
max_epochs_scvi = np.min([round((20000 / andata.n_obs) * 100), 100])
model_scvi.train(max_epochs=27)
andata.obsm["X_scVI"] = model_scvi.get_latent_representation()

## cluster data leiden
#overwrite the labels for one dataset with “unlabelled”
andata.obs['cluster'] = andata.obs['cluster'].cat.add_categories(["unlabelled"])
andata.obs.loc[andata.obs['batch'] == "1", "cluster"] = "unlabelled"
andata.obs['cluster'] = andata.obs['cluster'].astype(str)
# Then, if needed, re-convert it to a categorical column
andata.obs['cluster'] = andata.obs['cluster'].astype('category')


model_scanvi = scvi.model.SCANVI.from_scvi_model(
    model_scvi, labels_key='cluster', unlabeled_category="unlabelled"
)

max_epochs_scanvi = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
model_scanvi.train(max_epochs=max_epochs_scanvi)

andata.obsm["X_scANVI_lieden"] = model_scanvi.get_latent_representation()

## cluster data spatialleiden
#overwrite the labels for one dataset with “unlabelled”
andata.obs['spatialleiden'] = andata.obs['spatialleiden'].cat.add_categories(["unlabelled"])
andata.obs.loc[andata.obs['batch'] == "1", "spatialleiden"] = "unlabelled"
# Convert all values in the 'cluster' column to strings first
andata.obs['spatialleiden'] = andata.obs['spatialleiden'].astype(str)
# Then, if needed, re-convert it to a categorical column
andata.obs['spatialleiden'] = andata.obs['spatialleiden'].astype('category')

model_scanvi = scvi.model.SCANVI.from_scvi_model(
    model_scvi, labels_key='spatialleiden', unlabeled_category="unlabelled"
)

max_epochs_scanvi = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
model_scanvi.train(max_epochs=max_epochs_scanvi)

andata.obsm["X_scANVI_spatialleiden"] = model_scanvi.get_latent_representation()

andata_save = andata.copy()
andata_save.write_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_spatialleiden_harmony_scVI_scANVI.h5ad"))

