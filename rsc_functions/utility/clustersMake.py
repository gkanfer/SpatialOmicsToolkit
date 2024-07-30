from rsc_functions.utility.applyqc import applyqc 
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import seaborn as sns
import os
import gzip
import numpy as np
import rapids_singlecell as rsc

class clustersMake(applyqc):
    def __init__(self, leiden_resolution=0.5, leiden_iterations=2, *args, **kwargs):
        self.leiden_resolution = leiden_resolution
        self.leiden_iterations = leiden_iterations
        super().__init__(*args, **kwargs)
        self.perform_analysis()

    def perform_analysis(self):
        rsc.pp.normalize_total(self.andata)
        rsc.pp.log1p(self.andata)
        rsc.pp.scale(self.andata)
        rsc.pp.highly_variable_genes(self.andata,flavor="pearson_residuals", layer="counts")
        rsc.pp.pca(self.andata, n_comps=self.n_comps,random_state=self.random_state)
        rsc.pp.neighbors(self.andata,n_pcs=self.n_comps,use_rep='X_pca',n_neighbors=self.n_neighbors)
        rsc.tl.leiden(self.andata,random_state=self.random_state,resolution=self.leiden_resolution,key_added='cluster')
        return 