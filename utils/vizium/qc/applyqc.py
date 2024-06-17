from utils.vizium.viziumHD import viziumHD

class applyqc(viziumHD):
    '''
    A class to handle the processing and quality control reporting of Visium high-definition spatial transcriptomics data.
    
    Parameters
    ----------
    min_counts : int, optional
        Minimum counts for filtering cells. Default is 50.
    min_genes : int, optional
        Minimum genes for filtering cells. Default is 80.
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
    qcFilePrefix : str, optional
        Prefix for the QC report file. Default is an empty string.
    qc : bool, optional
        Flag to indicate whether QC should be performed. Default is True.
    '''
    def __init__(self,min_counts=50,min_genes = 80,rmvCellth = 50, mitoTh = 20, totalThr=10000, bins_total=40, bins_gene_plot=60, geneThr=4000,
                 bins_gene_plot_thr=60, qcFilePrefix="" ,mitoPlotPrefix = "",pcaPlotPrefix = "" ,n_comps = 50,qc=True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.min_counts = min_counts
        self.min_genes = min_genes
        self.rmvCellth = rmvCellth
        self.mitoTh = mitoTh
        self.totalThr = totalThr
        self.bins_total = bins_total
        self.bins_gene_plot = bins_gene_plot
        self.geneThr = geneThr
        self.bins_gene_plot_thr = bins_gene_plot_thr
        self.qcFilePrefix = qcFilePrefix
        self.mitoPlotPrefix = mitoPlotPrefix
        self.pcaPlotPrefix = pcaPlotPrefix
        self.n_comps = n_comps
        self.qc = qc
        if self.qc:
            self.apply_qc()
        
        
    def calcQCmat(self):
        '''
        Calculate QC metrics for the AnnData object.
        
        Returns
        -------
        AnnData
            The AnnData object with calculated QC metrics.
        '''
        self.andata.var_names_make_unique()
        self.andata.var["mt"] = self.andata.var_names.str.startswith("mt-")
        self.andata.var["ribo"] = self.andata.var_names.str.startswith(("RPS", "RPL"))
        self.andata.var["hb"] = self.andata.var_names.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(self.andata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
        return
        
    def apply_qc(self):
        '''
        Apply quality control to the AnnData object.
        '''
        self.andata = self.calcQCmat(self.andata)
        sc.pp.filter_cells(self.andata, min_counts=self.min_counts)
        sc.pp.filter_cells(self.andata, min_genes=self.min_genes)
        self.andata = self.andata[:, self.andata.var.n_cells_by_counts > self.rmvCellth]
        self.andata = self.andata[self.andata.obs["pct_counts_mt"] < self.mitoTh]
        self.qc_report()

    def qc_report(self):
        '''
        Generate a QC report and save it as a PDF.
        '''
        sc.pp.calculate_qc_metrics(self.andata, inplace=True)
        print("Starting QC report")
        with PdfPages(os.path.join(self.outPath, f'Quality_Control_{self.qcFilePrefix}.pdf')) as pdf:
            plt.rcParams['figure.dpi'] = 150
            plt.rcParams['font.family'] = ['serif']
            plt.rcParams['font.size'] = 12
            plt.rcParams['axes.labelsize'] = 12
            plt.rcParams['axes.titlesize'] = 12
            plt.rcParams['xtick.labelsize'] = 12
            plt.rcParams['ytick.labelsize'] = 12
            fig, axs = plt.subplots(2, 2, figsize=(8, 8))
            axs = axs.ravel()
            # Plot for total counts
            sns.histplot(self.andata.obs["total_counts"], kde=False, ax=axs[0])
            axs[0].set_title('Total Counts per Cell')
            axs[0].set_xlabel('Total Counts')
            axs[0].set_ylabel('Frequency')

            # Plot for total counts with a threshold
            sns.histplot(
                self.andata.obs["total_counts"][self.andata.obs["total_counts"] < self.totalThr],
                kde=False,
                bins=self.bins_total,
                ax=axs[1],
            )
            axs[1].set_title(f'Total Counts per Cell (Threshold < {self.totalThr})')
            axs[1].set_xlabel('Total Counts')
            axs[1].set_ylabel('Frequency')

            # Plot for number of genes by counts
            sns.histplot(self.andata.obs["n_genes_by_counts"], kde=False, bins=self.bins_gene_plot, ax=axs[2])
            axs[2].set_title('Number of Genes Detected per Cell')
            axs[2].set_xlabel('Number of Genes')
            axs[2].set_ylabel('Frequency')

            # Plot for number of genes by counts with a threshold
            sns.histplot(
                self.andata.obs["n_genes_by_counts"][self.andata.obs["n_genes_by_counts"] < self.geneThr],
                kde=False,
                bins=self.bins_gene_plot_thr,
                ax=axs[3],
            )
            axs[3].set_title(f'Number of Genes Detected per Cell (Threshold < {self.geneThr})')
            axs[3].set_xlabel('Number of Genes')
            axs[3].set_ylabel('Frequency')

            fig.tight_layout()
            pdf.savefig()
            plt.close()

    def plot_mitochondrial_data(self):
        '''
        Plot mitochondrial data for the AnnData object.
        '''
        with PdfPages(os.path.join(self.outPath, f'Quality_Control_mito{self.mitoPlotPrefix}.pdf')) as pdf:
            plt.rcParams['figure.dpi'] = 150
            plt.rcParams['font.family'] = ['serif']
            plt.rcParams['font.size'] = 12
            plt.rcParams['axes.labelsize'] = 12
            plt.rcParams['axes.titlesize'] = 12
            plt.rcParams['xtick.labelsize'] = 12
            plt.rcParams['ytick.labelsize'] = 12
            sc.pp.calculate_qc_metrics(self.andata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
            sns.jointplot(data=self.andata.obs, x="log1p_total_counts", y="log1p_n_genes_by_counts", kind="hex")
            pdf.savefig()
            plt.close()
            
    def perform_pca_and_save_plot(self):
        sc.tl.pca(self.andata, n_comps=n_comps, use_highly_variable=True)
        with PdfPages(os.path.join(self.outPath, f'Quality_Control_pca{self.pcaPlotPrefix }.pdf')) as pdf:
            plt.rcParams['figure.dpi'] = 150
            plt.rcParams['font.family'] = ['serif']
            plt.rcParams['font.size'] = 12
            plt.rcParams['axes.labelsize'] = 12
            plt.rcParams['axes.titlesize'] = 12
            plt.rcParams['xtick.labelsize'] = 12
            plt.rcParams['ytick.labelsize'] = 12
            plt.figure(figsize=(4, 3))
            n_pcs = 50

            # Plot the variance ratios and label the data points
            for i in range(n_pcs):
                x = i+1
                y = explained_variance_ratio[i]
                label = 'PC' + str(x)
                plt.plot(x, y, '')  # plots a blue circle marker
                plt.text(x, y, label, ha='right',fontsize = 4)  # positions text label at data point

            plt.xlabel('Principal Component')
            plt.ylabel('Variance Ratio')

            # Specify number of ticks on the y-axis
            ax = plt.gca()  # get current axes
            ax.yaxis.set_major_locator(ticker.MaxNLocator(10)) 
            
            pdf.savefig()
            plt.close()

            
            
            