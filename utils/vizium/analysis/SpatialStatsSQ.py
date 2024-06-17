from utils.vizium.analysis.analysis_sq import analysis_sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import seaborn as sns
import os
import gzip
import numpy as np
import scanpy as sc
import squidpy as sq

class SpatialStatsSQ(analysis_sq):
    def __init__(self,coord_type = "generic", *args, **kwargs):
        self.coord_type = coord_type
        super().__init__(*args, **kwargs)
        self.perform_spatial_statistics()
        
    def perform_spatial_statistics(self):
        '''
        Perform spatial statistics including neighborhood enrichment and co-occurrence.
        '''
        sq.gr.spatial_neighbors(self.andata, coord_type=self.coord_type, spatial_key="spatial")
        sq.gr.nhood_enrichment(self.andata, cluster_key="clusters")
        sq.gr.co_occurrence(self.andata, cluster_key="clusters")
        sq.gr.spatial_autocorr(self.andata, mode="moran")
        return

    