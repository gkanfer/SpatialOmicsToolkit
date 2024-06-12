import scanpy as sc
import rapids_singlecell as rsc
path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
adata = sc.read_visium(path=path_016)
rsc.get.anndata_to_GPU(adata)

