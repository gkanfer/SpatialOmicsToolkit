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
import voyagerpy as vp
from leidenalg import ModularityVertexPartition
import geopandas as gpd


class analysis_sq(applyqc):
    def __init__(self, leiden_resolution=0.5, leiden_iterations=2, *args, **kwargs):
        self.leiden_resolution = leiden_resolution
        self.leiden_iterations = leiden_iterations
        super().__init__(*args, **kwargs)
        self.perform_analysis()
        if self.voygerpyRead:
            self.perform_analysis_vp()
        else:
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
    
    def perform_analysis_vp(self):
        is_mt = self.andata.var_names.str.startswith('mt')
        vp.utils.add_per_cell_qcmetrics(self.andata, subsets={'mito': is_mt})

        spot_diameter_fullres = self.andata.uns['spatial']['scale'].pop('spot_diameter_fullres')
        self.andata.uns['spatial']['scale']['spot_diameter_fullres'] = {'pxl_col_in_fullres':spot_diameter_fullres,'pxl_row_in_fullres':spot_diameter_fullres}

        self.andata.uns['spatial']['scale']

        scale = 1
        scale_dict = self.andata.uns["spatial"].get("scale", {})
        spot_diam = scale_dict.get("spot_diameter_fullres")
        visium_spots = gpd.GeoSeries.from_xy(self.andata.obsm['spatial'][:,0], self.andata.obsm['spatial'][:,1]).scale(scale, scale, origin=(0, 0))

        self.andata.layers['counts'] = self.andata.X.copy()
        # Log-normalize the adata.X matrix
        vp.utils.log_norm_counts(self.andata, inplace=True)
        self.andata.layers['logcounts'] = self.andata.X.copy()

        self.andata.obs['var'] = np.var(self.andata.X.todense(), axis=1)

        gene_var = vp.utils.model_gene_var(self.andata.layers['logcounts'], gene_names=self.andata.var_names)
        hvgs = vp.utils.get_top_hvgs(gene_var)

        self.andata.var['highly_variable'] = False
        self.andata.var.loc[hvgs, 'highly_variable'] = True

        self.andata.X = vp.utils.scale(self.andata.X, center=True)
        sc.tl.pca(self.andata, use_highly_variable=True, n_comps=30, random_state=1337)
        self.andata.X = self.andata.layers['logcounts'].copy()

        self.andata.uns['config'] = OrderedDict()
        self.andata.uns["config"]["secondary_var_names"] = self.andata.var_names

        sc.pp.neighbors(
            self.andata,
            n_pcs=self.n_comps,
            use_rep='X_pca',
            method='gauss',
            n_neighbors=self.n_neighbors
        )
        sc.tl.leiden(
            self.andata,
            random_state=self.random_state,
            resolution=None,
            key_added='cluster',
            partition_type=ModularityVertexPartition
        ) 
        # de analysis
        
        markers = vp.utils.find_markers(self.andata, hvg = True)
        self.andata.var['symbol'] = self.andata.var['gene_ids'].values
        marker_genes = [
            marker.sort_values(by='p.value').iloc[0].name
            for _, marker in sorted(markers.items())]

        marker_genes_symbols = self.andata.var.loc[marker_genes, "symbol"].tolist()
        self.andata.var.loc[marker_genes, ["symbol"]]  
        self.andata.uns["marker_genes"] = marker_genes
        return
    

    



