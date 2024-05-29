from utils.vizium.viziumHD import viziumHD
import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import gzip
import numpy as np

class writeAnndataFile(viziumHD):
    def __init__(self,hdffileName,*args, **kwargs):
        self.hdffileName = hdffileName
        super().__init__(*args, **kwargs)
        
        
    def writeAnn(self):
        print("prepeare for writing")
        sc.pp.filter_cells(self.andata, min_counts = 50)
        sc.pp.filter_cells(self.andata, min_genes = 50)
        sc.pp.filter_cells(self.andata, max_counts = 1500)
        sc.pp.normalize_total(self.andata)
        sc.pp.log1p(self.andata)
        # sc.pp.scale(self.andata, max_value=10)
        sc.pp.pca(self.andata, n_comps = 10)
        # sc.pp.neighbors(self.andata)
        # sc.tl.umap(self.andata)
        # sc.tl.leiden(self.andata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)
        # # sq.pl.spatial_scatter(self.andata, color=["clusters"])
        self.andata.write(os.path.join(self.outPath, self.hdffileName))
        return self
        
    