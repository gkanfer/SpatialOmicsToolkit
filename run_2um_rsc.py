import scanpy as sc
#import squidpy as sq
import cupy as cp
import rapids_singlecell as rsc
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
import os
import gzip
import numpy as np

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)

path_002 = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_002um"
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out"


andata002 = sc.read_visium(path=path_002)
rsc.get.anndata_to_GPU(andata002)
rsc.pp.filter_cells(andata002, min_counts = 10)
andata002 = andata002[andata002.obs["pct_counts_mt"] < 20]
rsc.pp.normalize_total(andata002)
rsc.pp.log1p(andata002)
rsc.pp.scale(andata002, max_value=10)
andata002.obsm['spatial'] = np.array(andata002.obsm['spatial'], dtype=np.float64)

from matplotlib.backends.backend_pdf import PdfPages
# first preform PCA:
rsc.tl.pca(andata016, n_comps=50)
# then insted of using scanpy PCA ploting  
# rsc.pl.pca_variance_ratio(andata016, log=True, n_pcs=50)
# creat a custom PCA ploting

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Calculate the variance ratios
pca = andata002.uns['pca']
explained_variance_ratio = np.log(pca['variance_ratio'])


with PdfPages(os.path.join(pathout, f'Principal_Component_VR_2um.pdf')) as pdf:
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

rsc.pp.pca(andata002, n_comps=10)
rsc.pp.neighbors(andata002)
rsc.tl.umap(andata002)
rsc.tl.leiden(andata002, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)

andata002.write(os.path.join(pathout, "andata002_pca_neig_umap_leiden_file.h5ad"))

# with PdfPages(os.path.join(pathout, f'Single_Cell_Analysis_Conventional_2um_scatter.pdf')) as pdf:
#     # Set the plot parameters
#     plt.rcParams['figure.dpi'] = 150
#     plt.rcParams['font.family'] = ['serif']
#     plt.rcParams['font.size'] = 12
#     plt.rcParams['axes.labelsize'] = 12
#     plt.rcParams['axes.titlesize'] = 12
#     plt.rcParams['xtick.labelsize'] = 12
#     plt.rcParams['ytick.labelsize'] = 12

#     # Create the figure and axis
#     fig, ax = plt.subplots(1, 1, figsize=(4, 3))

#     # Plot the spatial scatter plot on the specified axis
#     sq.pl.spatial_scatter(andata002, color="clusters", ax=ax)

#     #     # Remove the legend
#     #     ax.get_legend().remove()

#     # Adjust layout and save to PDF
#     #     fig.tight_layout()
#     pdf.savefig()
#     plt.close()
