import scanpy as sc
import squidpy as sq
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
import os
import gzip
import numpy as np


def calcQCmat(andata):
    andata.var_names_make_unique()
    andata.var["mt"] = andata.var_names.str.startswith("mt-")
    andata.var["ribo"] = andata.var_names.str.startswith(("RPS", "RPL"))
    andata.var["hb"] = andata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(andata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    return andata


path_008 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_008um"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out"


andata008_ = sc.read_visium(path=path_008)
andata008 = calcQCmat(andata008_)
sc.pp.filter_cells(andata008, min_counts = 10)
andata008 = andata008[andata008.obs["pct_counts_mt"] < 20]
sc.pp.normalize_total(andata008)
sc.pp.log1p(andata008)
sc.pp.scale(andata008, max_value=10)
andata008.obsm['spatial'] = np.array(andata008.obsm['spatial'], dtype=np.float64)

from matplotlib.backends.backend_pdf import PdfPages
# first preform PCA:
sc.tl.pca(andata008, n_comps=50)
# then insted of using scanpy PCA ploting  
# sc.pl.pca_variance_ratio(andata016, log=True, n_pcs=50)
# creat a custom PCA ploting

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Calculate the variance ratios
pca = andata008.uns['pca']
explained_variance_ratio = np.log(pca['variance_ratio'])


with PdfPages(os.path.join(pathout, f'Principal_Component_VR_8um.pdf')) as pdf:
    # Create your own plot
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['font.family'] = ['serif']
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    
    fig, axs = plt.subplots(1, 1, figsize=(4, 3))

    n_pcs = 50

    # Plot the variance ratios and label the data points
    for i in range(n_pcs):
        x = i+1
        y = explained_variance_ratio[i]
        label = 'PC' + str(x)
        axs.plot(x, y, '')  # plots a blue circle marker
        axs.text(x, y, label, ha='right',fontsize = 4)  # positions text label at data point

    axs.set_xlabel('Principal Component')
    axs.set_ylabel('Variance Ratio')

    # Specify number of ticks on the y-axis
    ax = plt.gca()  # get current axes
    ax.yaxis.set_major_locator(ticker.MaxNLocator(10)) 
    fig.tight_layout()
    pdf.savefig()
    plt.close()

sc.pp.pca(andata008, n_comps=10)
sc.pp.neighbors(andata008)
sc.tl.umap(andata008)
sc.tl.leiden(andata008, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)

with PdfPages(os.path.join(pathout, f'Single_Cell_Analysis_Conventional_8um_scatter.pdf')) as pdf:
    # Set the plot parameters
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['font.family'] = ['serif']
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    # Create the figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Plot the spatial scatter plot on the specified axis
    sq.pl.spatial_scatter(andata008, color="clusters", ax=ax)

    #     # Remove the legend
    #     ax.get_legend().remove()

    # Adjust layout and save to PDF
    #     fig.tight_layout()
    pdf.savefig()
    plt.close()
