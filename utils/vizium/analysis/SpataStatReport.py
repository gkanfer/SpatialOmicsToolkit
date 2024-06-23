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


class SpataStatReport(SpatialStatsSQ):
    def __init__(self,n_genes_rank = 5 ,nbins_nhood = 100,sel_clusters = ['1','5'],reff_cluster = '1',num_genes_auto = 4, pv_cellphone = 0.001, alpha = 0.001, *args, **kwargs):
        self.n_genes_rank = n_genes_rank
        self.nbins_nhood = nbins_nhood 
        self.sel_clusters = sel_clusters
        self.reff_cluster = reff_cluster
        self.num_genes_auto = num_genes_auto
        self.pv_cellphone = pv_cellphone
        self.alpha = alpha
        super().__init__(*args, **kwargs)
        
    def report_rank_genes_groups(self):
        print("Starting Genes Ranking report")
        with PdfPages(os.path.join(self.outPath, f'Report_dotplotGeneGroups_{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            sc.pl.rank_genes_groups_dotplot(self.andata, groupby="clusters", standard_scale="var", n_genes=self.n_genes_rank, key="dea_clusters")
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            pdf.savefig()
            plt.close()
    
    def plot_nhood_enrichment(self):
        print("Starting neighborhood enrichment report")
        with PdfPages(os.path.join(self.outPath, f'Report__nhood_enr{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue -> White -> Red
            custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=self.nbins_nhood)
            fig, ax = plt.subplots()
            sq.pl.nhood_enrichment(self.andata, cluster_key="clusters", method="average", cmap=custom_cmap, vmin=-200, vmax=200, ax =ax)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            pdf.savefig()
            plt.close()
    
    def plot_spatial(self):
        print("Report spatial map")
        palette = sns.color_palette("tab20") + sns.color_palette("tab20b") + sns.color_palette("tab20c")
        listed_cmap = ListedColormap(palette)
        with PdfPages(os.path.join(self.outPath, f'Report_spatial_map_plot_{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            fig, ax = plt.subplots(1, 1, figsize=(4, 3))
            sq.pl.spatial_scatter(self.andata, color="clusters", img=False, ax=ax, palette=listed_cmap)
            fig.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            pdf.savefig()
            plt.close()
        
    
    def plot_spatial_clasterList(self):
        print("Report spatial gene map")
        palette = sns.color_palette("tab20") + sns.color_palette("tab20b") + sns.color_palette("tab20c")
        listed_cmap = ListedColormap([palette[clust] for clust in range(len(self.sel_clusters))])
        with PdfPages(os.path.join(self.outPath, f'Report_selectCluster_spatial_plot_{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            fig, ax = plt.subplots(1, 1, figsize=(4, 3))
            sq.pl.spatial_scatter(self.andata, groups=self.sel_clusters, color="clusters", img=False, ax=ax, palette=listed_cmap)
            fig.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            pdf.savefig()
            plt.close()
    
    def plot_co_occurrence(self):
        print("start coocurance analysis report")
        palette = sns.color_palette("tab20") + sns.color_palette("tab20b") + sns.color_palette("tab20c")
        listed_cmap = ListedColormap(palette)
        with PdfPages(os.path.join(self.outPath, f'Report__selectCluster_co_occurrence_{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            fig, ax = plt.subplots(1, 1, figsize=(4, 3))
            sq.pl.co_occurrence(self.andata, cluster_key="clusters", clusters=self.reff_cluster, palette=listed_cmap)
            plt.ylabel("co-occurrence ratio")
            fig.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            pdf.savefig()
            plt.close()
    
    def plot_top_autocorr(self):
        '''
        Plot top and bottom autocorrelation.
        '''
        # Ensure num_genes is within the allowed range
        print("start autocorlation report")
        num_genes = min(max(self.num_genes_auto, 1), 20)
        
        # Sort the genes based on Moran's I
        autocorr = (
            self.adata.uns["moranI"]["I"].sort_values(ascending=False).head(num_genes).index.tolist()
        )
        
        # Set up the plotting context
        sns.set_context("paper", font_scale=2)
        
        # Calculate the number of rows needed
        num_rows = (num_genes + 3) // 4  # Each row can have up to 4 plots
        
        with PdfPages(os.path.join(self.outPath, f'Report_Spatial_autocorr_{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            for i in range(num_rows):
                # Get the subset of genes to plot in the current row
                genes_subset = autocorr[i*4 : (i+1)*4]
            
            # Create the plot
            fig, axes = plt.subplots(1, len(genes_subset), figsize=(5 * len(genes_subset), 5))
            if len(genes_subset) == 1:
                axes = [axes]
            
            for ax, gene in zip(axes, genes_subset):
                sq.pl.spatial_scatter(
                    self.andata, color=gene, size=1, cmap="Reds", img=False, ax=ax, colorbar=False
                )
            
            plt.subplots_adjust(wspace=0.3)
            fig.tight_layout()
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
            # Save the current figure to the PDF
            pdf.savefig(fig)
            plt.close(fig)
            
    def run_cellphonedb(self):
        sq.gr.ligrec(self.andata,
        n_perms=1000,
        cluster_key="clusters",
        copy=False,
        use_raw=False,
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},)
        print("Starting Ligand Receptor interaction analysis")
        with PdfPages(os.path.join(self.outPath, f'Report__ligandReceptor{self.FilePrefix}.pdf')) as pdf:
            self.set_image_para()
            sq.pl.ligrec(self.andata, cluster_key="clusters",pvalue_threshold = 0.001,remove_empty_interactions = True, alpha = 0.001)
            