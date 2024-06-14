from skimage import io
import numpy as np
import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import os
import gzip
import numpy as np
import celltypist
from celltypist import models

colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]  # Blue -> White -> Red
n_bins = 100  # Discretize into 100 bins
custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=n_bins)

# Define the path to the saved .h5ad file
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out"
file_path = os.path.join(pathout, 'leiden0.5_rank_16um.h5ad')

# Load the .h5ad file
andata016 = sc.read(file_path)

sq.gr.spatial_neighbors(andata016, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(andata016, cluster_key="clusters")

with PdfPages(os.path.join(pathout, f'nhood_enrichment_generic_16um.pdf')) as pdf:
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['font.family'] = ['serif']
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    fig, ax = plt.subplots()  # Increase figure size

    sq.pl.nhood_enrichment(
        andata016,
        title = "Neighborhood Enrichment generic",
        cluster_key="clusters",
        method = "average",
        cmap = custom_cmap,
        vmin = -200,
        vmax = 200,
        ax = ax)

    fig.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)  # Adjust subplot parameters

    pdf.savefig(fig)
    plt.close()
    
sq.gr.spatial_neighbors(andata016,coord_type="grid", n_neighs=6, n_rings=8, key_added='spatial_neighbors')
sq.gr.nhood_enrichment(andata016, cluster_key="clusters")

with PdfPages(os.path.join(pathout, f'nhood_enrichment_grid_6_8_16um.pdf')) as pdf:
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['font.family'] = ['serif']
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    fig, ax = plt.subplots()  # Increase figure size

    sq.pl.nhood_enrichment(
        andata016,
        title = "Neighborhood Enrichment grid",
        cluster_key="clusters",
        method="average",
        cmap=custom_cmap,
        vmin=-200,
        vmax=200,
        ax = ax)

    fig.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)  # Adjust subplot parameters

    pdf.savefig(fig)
    plt.close()

