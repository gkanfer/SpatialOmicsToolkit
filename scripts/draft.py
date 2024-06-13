from skimage import io
import numpy as np
import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import gzip
import numpy as np
import celltypist
from celltypist import models

'''
This script performs quality control, normalization, scaling, PCA, UMAP, Leiden clustering, differential expression analysis, and visualizes the results in a dot plot for spatial omics data using Scanpy and Squidpy. The output is saved as a PDF and the processed data is saved in an .h5ad file.
'''

plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

# Define the path to the saved .h5ad file
pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out"
file_path = os.path.join(pathout, 'example_file.h5ad')

# Load the .h5ad file
andata016 = sc.read(file_path)

# Verify the loaded data
print(andata016)