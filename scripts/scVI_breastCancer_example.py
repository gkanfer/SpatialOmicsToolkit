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
andata = sc.read_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_harmony.h5ad"))

scvi.model.SCVI.setup_anndata(andata, layer="counts", batch_key="batch")
model_scvi = scvi.model.SCVI(andata)
max_epochs_scvi = np.min([round((20000 / andata.n_obs) * 100), 100])
model_scvi.train(max_epochs=27)
andata.obsm["X_scVI"] = model_scvi.get_latent_representation()
andata_save = adata_concat.copy()
#andata_save.write_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_harmony_scVI.h5ad"))




# andata.obsm["X_scVI"].shape
# (293960, 10)
#
# sc.pp.neighbors(andata, use_rep="X_scVI",key_added = 'scVI')
# sc.tl.leiden(andata, random_state=1337, resolution=0.5, key_added='cluster_scVI', neighbors_key='scVI')
# sc.tl.umap(adata_concat, neighbors_key="scVI")
# adata_concat.obsm["scVI_umap"] = adata_concat.obsm["X_umap"]
# andata_save = adata_concat.copy()
# andata_save.write_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_harmony_scVI.h5ad"))



