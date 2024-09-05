import cupy as cp
import cupyx
import scanpy as sc
import scanpy as sc
import numpy as np
from cupyx.scipy.sparse import csr_matrix
import os
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import random
import pandas as pd
import voyagerpy as vp
import geopandas as gpd
import libpysal as lps
from collections import OrderedDict
import scipy.sparse as sp
from scipy.sparse import csr_matrix
def run(sample_size):
    path_016 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Human_Colon_Cancer_binned_outputs/binned_outputs/square_016um"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
    andata = sc.read_h5ad(os.path.join(pathout, "andata_save_colon.h5ad"))
    andata.uns['config'] = OrderedDict()
    andata.uns['config'] = OrderedDict()
    andata.uns["config"]["secondary_var_names"] = andata.var_names

    import os
    import pickle

    file_path = os.path.join(pathout, "andata_uns_mtracies__colon.pkl")

    if os.path.getsize(file_path) > 0:
        with open(file_path, 'rb') as buff:
            andata_uns_mtracies = pickle.load(buff)
    else:
        print("File is empty. Cannot load data.")
        andata_uns_mtracies = None
    
    andata_sub = andata.copy()
    andata_sub.X = csr_matrix(andata_sub.X)
    andata_sub = sc.pp.subsample(andata_sub, n_obs=sample_size,copy=True) 
    
    sparse_dist_matrix = andata_sub.obsp['distances'].tocsr()
    sparse_inv_matrix = sparse_dist_matrix.copy()
    sparse_inv_matrix.data = 1 / sparse_inv_matrix.data
    sparse_inv_matrix.data[sparse_inv_matrix.data == float('inf')] = 0
    

sample_size = [5000,10_000,20_000,50_000,60_000]
