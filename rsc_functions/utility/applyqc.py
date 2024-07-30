from rsc_functions.StDatareader_rsc import StDatareader_rsc 
from collections import OrderedDict
import numpy as np
import rapids_singlecell as rsc


class applyqc(StDatareader_rsc):
    '''
    A class to handle the processing and quality control reporting of Visium high-definition spatial transcriptomics data.
    
    Parameters
    ----------
    min_counts : int
        Minimum counts for filtering cells. Default is 50. 0-inf
    min_genes : int
        Minimum genes for filtering cells. Default is 80. 0-inf
    min_gene_count : int
        Minimum genes for filtering genes. Default is 0. 0-inf
    max_gene_count : int
        Maximum genes for filtering genes. Default is 0. 0-inf
    rmvCellth : int, optional
        Threshold for removing cells based on the number of cells by counts. Default is 50.
    mitoTh : int, optional
        Threshold for the percentage of mitochondrial counts. Default is 20.
    totalThr : int, optional
        Threshold for filtering out cells based on total counts below this threshold. Default is 10000.
    bins_total : int, optional
        Number of bins to use for the histogram of total counts. Default is 40.
    bins_gene_plot : int, optional
        Number of bins to use for the histogram plotting the number of genes by counts. Default is 60.
    geneThr : int, optional
        Threshold for filtering out cells based on the number of genes detected below this threshold. Default is 4000.
    bins_gene_plot_thr : int, optional
        Number of bins to use for the histogram of the number of genes by counts with a threshold. Default is 60.
    voygerpyRead : bool, optional
        Whether to read data using voygerpy. Default is False.
    dataset_key : str, optional
        Key or identifier for the dataset. Default is "Visium_HD_Mouse_Brain".
    '''
    def __init__(self,mitochon = "mt-", min_counts=50,min_genes = 80,min_gene_count = 0,max_gene_count = 0, rmvCellth = 50, mitoTh = 20, totalThr = 10000, bins_total = 40, bins_gene_plot = 60, geneThr = 4000, bins_gene_plot_thr = 60, n_comps = 15, voygerpyRead = False, dataset_key = "Visium_HD_Mouse_Brain", n_neighbors = 30, random_state = 34566, *args, **kwargs):
        self.mitochon = mitochon #Mitochondia-encoded genes (gene names start with prefix mt- or MT-)
        self.min_counts = min_counts
        self.min_genes = min_genes
        self.min_gene_count = min_gene_count
        self.max_gene_count = max_gene_count
        self.rmvCellth = rmvCellth
        self.mitoTh = mitoTh
        self.totalThr = totalThr
        self.bins_total = bins_total
        self.bins_gene_plot = bins_gene_plot
        self.geneThr = geneThr
        self.bins_gene_plot_thr = bins_gene_plot_thr
        self.n_comps = n_comps
        self.voygerpyRead = voygerpyRead
        self.dataset_key = dataset_key
        self.n_neighbors = n_neighbors
        self.random_state = random_state
        super().__init__(*args, **kwargs)
        self.apply_qc()    
        
    
    def apply_qc(self):
        '''
        Apply quality control to the AnnData object.
        '''
        self.andata.var_names_make_unique()
        self.andata.obsm['spatial'] = np.array(self.andata.obsm['spatial'], dtype=np.float64)
        self.andata.uns['config'] = OrderedDict()
        self.andata.uns["config"]["secondary_var_names"] = self.andata.var_names
        rsc.pp.flag_gene_family(self.andata, gene_family_name="MT", gene_family_prefix="mt-")
        rsc.pp.calculate_qc_metrics(self.andata, qc_vars=["MT"])
        if self.min_counts > 0:
            rsc.pp.filter_cells(self.andata, min_count=self.min_counts,qc_var = 'total_counts')
        if self.min_gene_count > 0:
            rsc.pp.filter_genes(self.andata, min_count=self.min_gene_count)
        if self.max_gene_count > 0:
            rsc.pp.filter_genes(self.andata, max_count=self.max_gene_count)
        self.andata.layers['counts'] = self.andata.X.copy()
        return

    
            