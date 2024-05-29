from utils.vizium.viziumHD import viziumHD
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import gzip
import numpy as np


class addQCplots(viziumHD):
    def __init__(self, qc_add,qcFilePrefixAddQCreport, marker="Mt-", n_top_genes=2000, *args, **kwargs):
        '''
        This class extends viziumHD and provides methods to perform quality control analysis and generate additional QC plots.
        
        Parameters
        ----------
        qc_add : list
        List of additional QC plot functions to generate.
        marker : str, optional
        Prefix for genes to identify specific gene types. Default is "Mt-".
        n_top_genes : int, optional
        Number of top genes to consider in the analysis. Default is 2000.
        '''
        self.qc_add = qc_add
        self.marker = marker
        self.n_top_genes = n_top_genes
        self.qcFilePrefixAddQCreport = self.qcFilePrefixAddQCreport
        super().__init__(*args, **kwargs)

    def qc_metrics(self):
        self.andata.var["mt"] = self.andata.var_names.str.startswith(self.marker)
        self.andata.var["ribo"] = self.andata.var_names.str.startswith(("RPS", "RPL"))
        self.andata.var["hb"] = self.andata.var_names.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(
            self.andata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
        )
        sns.jointplot(data=self.andata.obs,x="log1p_total_counts",y="log1p_n_genes_by_counts",kind="hex",)

    def rm_high_variable_genes(self):
        sc.pp.highly_variable_genes(self.andata, n_top_genes=self.n_top_genes)
        sc.pl.highly_variable_genes(self.andata)

    def quality_control_pca(self):
        sc.pp.pca(self.andata)
        sc.pl.pca_variance_ratio(self.andata, n_pcs=50, log=True)
    
    def add_additional_qc_plots(self):
        qc_file = os.path.join(self.outPath, f'{self.qcFilePrefix}.pdf')
        if not os.path.exists(qc_file):
            self.qcReport()

        with PdfPages(qc_file) as pdf:
            plt.rcParams['figure.dpi'] = 150
            plt.rcParams['font.family'] = ['serif']
            plt.rcParams['font.size'] = 12
            plt.rcParams['axes.labelsize'] = 12
            plt.rcParams['axes.titlesize'] = 12
            plt.rcParams['xtick.labelsize'] = 12
            plt.rcParams['ytick.labelsize'] = 12

            for plot_func in self.qc_add:
                if hasattr(self, plot_func):
                    fig = plt.figure()
                    getattr(self, plot_func)()
                    pdf.savefig(fig)
                    plt.close(fig)