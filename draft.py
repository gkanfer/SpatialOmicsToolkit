from spatialdata_io import xenium
import spatialdata as sd

xenium_path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/outs"
path_write = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/zarr_rep1"

sdata = xenium(
    path=str(xenium_path),
    n_jobs=8,
    cells_boundaries=True,
    nucleus_boundaries=True,
    morphology_focus=True,
    cells_as_circles=False,
)
print("done")