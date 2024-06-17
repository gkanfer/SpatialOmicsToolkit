from utils.vizium.qc.applyqc import applyqc

class analysis_qc(applyqc):
    def __init__(self, n_comps=10, leiden_resolution=0.5, leiden_iterations=2,min_in_group_fraction = 0.01,max_out_group_fraction = 0.01 ,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_comps = n_comps
        self.leiden_resolution = leiden_resolution
        self.leiden_iterations = leiden_iterations
        self.perform_analysis()
        
    def perform_analysis(self):
        '''
        Perform data normalization, transformation, PCA, UMAP, and clustering.
        '''
        sc.pp.normalize_total(self.andata)
        sc.pp.log1p(self.andata)
        sc.pp.highly_variable_genes(self.andata)
        sc.pp.scale(self.andata)
        sc.pp.pca(self.andata, n_comps=self.n_comps)
        sc.pp.neighbors(self.andata)
        sc.tl.umap(self.andata)
        sc.tl.leiden(self.andata, key_added='clusters', flavor="igraph", directed=False, resolution=self.leiden_resolution, n_iterations=self.leiden_iterations)
        sc.tl.rank_genes_groups(self.andata, groupby="clusters", method="wilcoxon", key_added="dea_clusters")
        
    def report_rank_genes_groups(self):
        with PdfPages(os.path.join(self.outPath, f'Report_dotplotGeneGroups_{self.PlotPrefix}.pdf')) as pdf:
            sc.pl.rank_genes_groups_dotplot(self.andata, groupby="clusters", standard_scale="var", n_genes=5, key="dea_clusters")
            pdf.savefig()
            plt.close()

    



