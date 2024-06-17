from utils.vizium.analysis.analysis_sq import analysis_sq

class SpatialStatsSQ(analysis_sq):
    def __init__(self, custom_cmap=None, cluster_to_plot='5', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.perform_spatial_statistics()
        
    def perform_spatial_statistics(self):
        '''
        Perform spatial statistics including neighborhood enrichment and co-occurrence.
        '''
        sq.gr.spatial_neighbors(self.andata, coord_type="generic", spatial_key="spatial")
        sq.gr.nhood_enrichment(self.andata, cluster_key="clusters")
        sq.gr.co_occurrence(self.andata, cluster_key="clusters")
        sq.gr.spatial_autocorr(self.andata, mode="moran")
        
class SpataStatReport(SpatialStatsSQ):
    def __init__(self,nbins_nhood = 100,sel_clusters = ['1','5'],reff_cluster = '1',num_genes_auto = 4):
        self.nbins_nhood = nbins_nhood 
        self.sel_clusters = sel_clusters
        self.reff_cluster = reff_cluster
        self.num_genes_auto = num_genes_auto
        
    def plot_nhood_enrichment(self):
        '''
        Plot neighborhood enrichment.
        
        Parameters
        ----------
        save_path : str, optional
            Path to save the neighborhood enrichment plot. Default is "nhood_enrichment.pdf".
        '''
        with PdfPages(os.path.join(self.outPath, f'Report__nhood_enr{self.PlotPrefix}.pdf')) as pdf:
            colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue -> White -> Red
            custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=self.nbins_nhood)
            sq.pl.nhood_enrichment(self.andata, cluster_key="clusters", method="average", cmap=custom_cmap, vmin=-200, vmax=200, figsize=(3, 3))
            pdf.savefig()
            plt.close()
    
    def plot_spatial_clasterList(self):
        listed_cmap = ListedColormap([clast for clast in self.sel_clusters])
        with PdfPages(os.path.join(self.outPath, f'Report_selectCluster_spatial_plot_{self.PlotPrefix}.pdf')) as pdf:
            fig, ax = plt.subplots(1, 1, figsize=(4, 3))
            sq.pl.spatial_scatter(andata016, groups=sel_cluters, color="clusters",img=False, ax=ax, palette=listed_cmap)
            pdf.savefig()
            plt.close()
    
    
    
    def plot_co_occurrence(self, save_path="co_occurrence.pdf"):
        '''
        Plot co-occurrence.
        
        Parameters
        ----------
        save_path : str, optional
            Path to save the co-occurrence plot. Default is "co_occurrence.pdf".
        '''
        palette = sns.color_palette("tab20") + sns.color_palette("tab20b") + sns.color_palette("tab20c")
        listed_cmap = ListedColormap(palette)
        with PdfPages(os.path.join(self.outPath, f'Report__selectCluster_co_occurrence_{self.PlotPrefix}.pdf')) as pdf:
            sq.pl.co_occurrence(self.andata, cluster_key="clusters", clusters=self.reff_cluster, figsize=(4, 3), palette=listed_cmap)
            plt.ylabel("co-occurrence ratio")
            pdf.savefig()
            plt.close()
    
    def plot_top_bottom_autocorr(self, num_view=4, save_path="spatial_autocorr.pdf"):
        '''
        Plot top and bottom autocorrelation.
        
        Parameters
        ----------
        num_view : int, optional
            Number of top and bottom autocorrelated genes to view. Default is 4.
        save_path : str, optional
            Path to save the autocorrelation plot. Default is "spatial_autocorr.pdf".
        '''
        # Ensure num_genes is within the allowed range
        num_genes = min(max(self.num_genes_auto, 1), 20)
        
        # Sort the genes based on Moran's I
        top_autocorr = (
            adata.uns["moranI"]["I"].sort_values(ascending=False).head(num_genes).index.tolist()
        )
        
        # Set up the plotting context
        sns.set_context("paper", font_scale=2)
        
        # Calculate the number of rows needed
        num_rows = (num_genes + 3) // 4  # Each row can have up to 4 plots
        
        with PdfPages(os.path.join(self.outPath, f'Report_Spatial_autocorr_{self.PlotPrefix}.pdf')) as pdf:
            for i in range(num_rows):
            # Get the subset of genes to plot in the current row
            genes_subset = top_autocorr[i*4 : (i+1)*4]
            
            # Create the plot
            fig, axes = plt.subplots(1, len(genes_subset), figsize=(5 * len(genes_subset), 5))
            if len(genes_subset) == 1:
                axes = [axes]
            
            for ax, gene in zip(axes, genes_subset):
                sq.pl.spatial_scatter(
                    adata, color=gene, size=1, cmap="Reds", img=False, ax=ax, colorbar=False
                )
            
            plt.subplots_adjust(wspace=0.3)
            
            # Save the current figure to the PDF
            pdf.savefig(fig)
            plt.close(fig)
