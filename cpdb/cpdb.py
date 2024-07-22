import numpy as np
import geopandas as gpd
import pandas as pd

import scanpy as sc
import squidpy as sq
import voyagerpy as vp
import seaborn as sns
import os
import pickle
from matplotlib.pyplot import imread
from collections import OrderedDict
import json
from matplotlib import pyplot as plt
    
from cellphonedb.src.core.methods import cpdb_analysis_method
    
path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
# andata016_ = sc.read_visium(path=path_016)
# andata016_
# andata016_.var_names_make_unique()
# andata016_.obsm['spatial'] = np.array(andata016_.obsm['spatial'], dtype=np.float64)


# sc.pp.filter_cells(andata016_, min_counts=1000)

# sc.pp.filter_cells(andata016_, min_genes=1000)

# sc.pp.filter_genes(andata016_, min_counts=1000)

# sc.pp.filter_genes(andata016_, max_counts=6281)

# andata016_.obsm['spatial'] = np.array(andata016_.obsm['spatial'], dtype=np.float64)
# andata016_.uns['spatial']['img'] = andata016_.uns['spatial']['Visium_HD_Mouse_Brain'].pop("images")
# andata016_.uns['spatial']['scale'] = andata016_.uns['spatial']['Visium_HD_Mouse_Brain'].pop("scalefactors")
# andata016_.uns['spatial']['metadata'] = andata016_.uns['spatial']['Visium_HD_Mouse_Brain'].pop("metadata")
# andata016_.uns['spatial'].pop("Visium_HD_Mouse_Brain")

# # change order of images
# images = andata016_.uns['spatial'].pop('img')
# images_hires = {'lowres':images['lowres'],'hires':images['hires']}
# andata016_.uns['spatial']['img'] = images_hires


# is_mt = andata016_.var_names.str.startswith('mt')
# vp.utils.add_per_cell_qcmetrics(andata016_, subsets={'mito': is_mt})


# spot_diameter_fullres = andata016_.uns['spatial']['scale'].pop('spot_diameter_fullres')
# andata016_.uns['spatial']['scale']['spot_diameter_fullres'] = {'pxl_col_in_fullres':spot_diameter_fullres,'pxl_row_in_fullres':spot_diameter_fullres}
# # insted of vp.spatial.get_visium_spots(andata016_, with_radius=False) I have done:
# #scale = andata016_.uns['spatial']['scale']['tissue_lowres_scalef']
# scale = 1
# scale_dict = andata016_.uns["spatial"].get("scale", {})
# spot_diam = scale_dict.get("spot_diameter_fullres")
# visium_spots = gpd.GeoSeries.from_xy(andata016_.obsm['spatial'][:,0], andata016_.obsm['spatial'][:,1]).scale(scale, scale, origin=(0, 0))
# _ = vp.spatial.set_geometry(andata016_, geom="spot_poly", values=visium_spots)
# andata016_.uns['config'] = OrderedDict()
# andata016_.uns["config"]["secondary_var_names"] = andata016_.var_names
# pathout = "/data/kanferg/Sptial_Omics/VoyagerPy_fork/voyagerpy/out"

# qc_features = ["sum", "detected", "subsets_mito_percent"]
# andata016_.uns['config'] = OrderedDict()
# andata016_.uns["config"]["secondary_var_names"] = andata016_.var_names   

# # The original count data
# andata016_.layers['counts'] = andata016_.X.copy()
# # Log-normalize the adata.X matrix
# vp.utils.log_norm_counts(andata016_, inplace=True)
# andata016_.layers['logcounts'] = andata016_.X.copy()


# gene_var = vp.utils.model_gene_var(andata016_.layers['logcounts'], gene_names=andata016_.var_names)
# hvgs = vp.utils.get_top_hvgs(gene_var)

# andata016_.var['highly_variable'] = False
# andata016_.var.loc[hvgs, 'highly_variable'] = True

# andata016_.X = vp.utils.scale(andata016_.X, center=True)
# sc.tl.pca(andata016_, use_highly_variable=True, n_comps=30, random_state=1337)
# andata016_.X = andata016_.layers['logcounts'].copy()

# from leidenalg import ModularityVertexPartition
# sc.pp.neighbors(
#     andata016_,
#     n_pcs=9,
#     use_rep='X_pca',
#     method='gauss',
#     n_neighbors=80
# )
# sc.tl.leiden(
#     andata016_,
#     random_state=29,
#     resolution=None,
#     key_added='cluster',
#     partition_type=ModularityVertexPartition
# )


cpdb_file_path = "/data/kanferg/cellphonedb/NatureProtocols2024_case_studies/v5.0.0/cellphonedb.zip"
meta_file_path = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2/test_meta.pickle"
counts_file_path = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2/test_counts.pickle"
out_path = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/cpdb_out" 

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(cpdb_file_path = cpdb_file_path, meta_file_path = meta_file_path, counts_file_path = counts_file_path, counts_data = 'ensembl', output_path = out_path) 

 