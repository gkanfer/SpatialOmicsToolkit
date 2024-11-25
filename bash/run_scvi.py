import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scvi
import torch
from scvi.autotune import ModelTuner
from ray import tune
import ray

torch.set_float32_matmul_precision("high")

print("CUDA available:", torch.cuda.is_available())

pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_4"
andata_combined = sc.read_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_harmony_scVI_scANVI_unintigrated.h5ad"))
andata_bc = andata_combined[andata_combined.obs['sample']=='ST'].copy()

model_cls = scvi.model.SCVI
model_cls.setup_anndata(andata_bc, categorical_covariate_keys = ['batch'], continuous_covariate_keys=['percent_mito'])
tuner = ModelTuner(model_cls)

model = model_cls(andata_bc)

search_space = {
    "n_hidden": tune.choice([92, 128]),
    "n_latent": tune.choice([10, 20, 30, 40, 50, 60]),
    #"n_layers": tune.choice([1, 2, 3]),
    "lr": tune.loguniform(1e-4, 1e-2),
    "gene_likelihood": tune.choice(["nb", "zinb"])
}

# Specify a storage path (e.g., a local directory for Ray's outputs)
#run_config = RunConfig(storage_path="./ray_results")

# Run the tuner with the updated configuration
results = tuner.fit(
    andata_bc,
    metric="validation_loss",
    resources={'gpu': 3},  # specify GPU resources
    search_space=search_space,
    num_samples=10,
    max_epochs=2,)

best_vl = 10000
best_i = 0
for i, res in enumerate(results.results):
    vl = res.metrics['validation_loss']

    if vl < best_vl:
        best_vl = vl
        best_i = i
        
results.results[best_i]

print(f'{results.results[best_i]}')