from spatialomicstoolkit.analysis.analysis_sq import analysis_sq
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
from scipy.sparse import csr_matrix
import voyagerpy as vp

class SpatialStatsSQ(analysis_sq):
    '''
    wn is wight Matrix 
    '''
    def __init__(self,coord_type = "generic", *args, **kwargs):
        self.coord_type = coord_type
        super().__init__(*args, **kwargs)
        if self.voygerpyRead:
            self.perform_spatial_statistics_vp()
        else:
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

    def perform_spatial_statistics_vp(self):
        sc.pp.neighbors(
            self.andata,
            n_neighbors=self.n_neighbors,
            n_pcs=self.n_comps,
            use_rep='X_pca',
            knn=True,
            random_state=self.random_state,
            method='gauss', # one of umap, gauss, rapids
            metric='euclidean', # many metrics available,
            key_added='knn')
        dist = self.andata.obsp['knn_distances'].copy()
        dist.data = 1 / dist.data

        # row normalize the matrix, this makes the matrix dense.
        dist /= dist.sum(axis=1)

        # convert dist back to sparse matrix
        self.andata.obsp["knn_weights"] = csr_matrix(dist)

        del dist
        knn_graph = "knn_weights"
        
        self.andata.obsp["knn_connectivities"] = (self.andata.obsp[knn_graph] > 0).astype(int)
        vp.spatial.set_default_graph(self.andata, knn_graph)
        vp.spatial.to_spatial_weights(self.andata, graph_name=knn_graph)

        qc_features = ["sum", "detected", "subsets_mito_percent"]
        morans = vp.spatial.moran(self.andata, qc_features, graph_name=knn_graph)
        self.andata.uns['spatial']['moran'][knn_graph].loc[qc_features, ["I"]]


        vp.spatial.compute_spatial_lag(
            self.andata,
            qc_features,
            graph_name=knn_graph,
            inplace=True
        )
        vp.spatial.local_moran(self.andata, qc_features, graph_name=knn_graph)
        return



    