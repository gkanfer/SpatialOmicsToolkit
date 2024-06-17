from utils.vizium.qc.applyqc import applyqc
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import seaborn as sns
import os
import gzip
import numpy as np


class analysis_sq(applyqc):
    def __init__(self, leiden_resolution=0.5, leiden_iterations=2, *args, **kwargs):
        self.leiden_resolution = leiden_resolution
        self.leiden_iterations = leiden_iterations
        super().__init__(*args, **kwargs)
        self.perform_analysis()
        
    def perform_analysis(self):
        '''
        Perform data normalization, transformation, PCA, UMAP, and clustering.
        '''
        print("intiate processing")
        sc.pp.normalize_total(self.andata)
        sc.pp.log1p(self.andata)
        sc.pp.highly_variable_genes(self.andata)
        sc.pp.scale(self.andata)
        print("intiate spatial analysis")
        sc.pp.pca(self.andata, n_comps=self.n_comps)
        sc.pp.neighbors(self.andata)
        sc.tl.umap(self.andata)
        sc.tl.leiden(self.andata, key_added='clusters', flavor="igraph", directed=False, resolution=self.leiden_resolution, n_iterations=self.leiden_iterations)
        sc.tl.rank_genes_groups(self.andata, groupby="clusters", method="wilcoxon", key_added="dea_clusters")
        return 
        
    

    



