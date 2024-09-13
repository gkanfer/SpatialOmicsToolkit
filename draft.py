from numpy.random import default_rng

import matplotlib.pyplot as plt

import scanpy as sc
import squidpy as sq
from anndata import AnnData


import os
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import random
import pandas as pd


pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/R_simulation"
counts = pd.read_csv(os.path.join(pathout,'Rsim_count.csv'))
adata = AnnData(counts.iloc[:,1:].to_numpy(), obsm={"spatial": pd.read_csv(os.path.join(pathout,'Rsim_centroids.csv')).loc[:,['x','y']].to_numpy()})
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=2)
sc.pp.neighbors(adata,n_neighbors = 10,method='gauss')
