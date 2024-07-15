from utils.vizium.analysis.SpatialStatsSQ import SpatialStatsSQ
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import os
import gzip
import numpy as np
import scanpy as sc
import squidpy as sq

class LigandReceptorAna(SpatialStatsSQ):
    def __init__(self,pv_cellphone = 0.001,alpha = 0.001,*args, **kwargs):
        self.pv_cellphone = pv_cellphone
        self.alpha = alpha
        super().__init__(*args, **kwargs)
        
    def run_cellphonedb(self):
        sq.gr.ligrec(self.andata,
        n_perms=1000,
        cluster_key="cluster",
        copy=False,
        use_raw=False,
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},)
        print("Starting Ligand Receptor interaction analysis")
        with PdfPages(os.path.join(self.outPath, f'Report__ligandReceptor{self.FilePrefix}.pdf')) as pdf:
            sq.pl.ligrec(self.andata, cluster_key="cluster",pvalue_threshold = 0.1,remove_empty_interactions = True, alpha = 0.01)
            