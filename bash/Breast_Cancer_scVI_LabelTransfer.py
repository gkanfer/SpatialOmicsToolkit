# This code uses the scvitools/1.2.0.gpu module on the Biowulf cluster with 3 GPUs (A100).
# On average, one training epoch takes 30 seconds when using this module.
# However, it does not perform well with the version of scVI installed manually.
# With the manually installed scVI, one training epoch takes only 0.5 seconds.

# The tuning process is recommended for 100 samples (`num_samples`), 
# which involves training 100 different models, each with varying parameters.
# If tuning is done at maximum capacity (100 models, 100 epochs each), 
# the total time required would be:
# (0.5 seconds/epoch * 100 epochs/model * 100 models) / 60 = 83 hours.

# This duration is too long to be practical. 
# It is advisable to reduce the number of samples to 20, which would result in:
# (0.5 seconds/epoch * 100 epochs/model * 20 models) / 60 = 16 hours.
# This is a more reasonable runtime for tuning.

# For more details on the tuning process, refer to the official scVI documentation:
# https://docs.scvi-tools.org/en/1.2.0/tutorials/notebooks/tuning/autotune_scvi.html

import ray
import scanpy as sc
import scvi
import seaborn as sns
import torch
from ray import tune
from scvi import autotune
import os

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_4"
andata_combined = sc.read_h5ad(os.path.join(pathout, "adata_concat_BreastCancer_harmony_scVI_scANVI_unintigrated.h5ad"))
andata_bc = andata_combined[andata_combined.obs['sample']=='ST'].copy()

model_cls = scvi.model.SCVI
model_cls.setup_anndata(andata_bc, categorical_covariate_keys = ['batch'])


search_space = {
    "model_params": {"n_hidden": tune.choice([92, 128]),
        "n_latent": tune.choice([10, 20, 30, 40, 50, 60]),
        "n_layers": tune.choice([1, 2, 3]),
        #"lr": tune.loguniform(1e-4, 1e-2),
        "gene_likelihood": tune.choice(["nb", "zinb"])},
    "train_params": {"max_epochs": 100},
}

ray.init(log_to_driver=False)
results = autotune.run_autotune(
    model_cls,
    data=andata_bc,
    mode="min",
    metrics="validation_loss",
    search_space=search_space,
    num_samples=20,
    resources={"gpu": 3})

df = results.result_grid.get_dataframe(filter_metric="accuracy", filter_mode="max")
df.to_csv(os.path.join(pathout,"tune_hyper_scVI.csv"))