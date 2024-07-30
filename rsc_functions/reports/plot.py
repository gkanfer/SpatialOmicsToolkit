import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import seaborn as sns
import os
import gzip
import numpy as np
import rapids_singlecell as rsc


def set_image_para(self):
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['font.family'] = ['serif']
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12



    
def plot_dist(andata,column,ax,type = 'obs', bins = 'auto',title = '',xlab = '',ylab =''):
    '''
    You can replace 'auto' with any other method (e.g., 'fd', 'doane', 'scott', 'rice', 'sturges', or 'sqrt')
    '''
    palette1 = sns.color_palette("colorblind",9)
    if type == 'obs':
        arr = andata.obs['column'].values
    else:
        arr = andata.var['column'].values
    bin_edges = np.histogram_bin_edges(arr, bins='auto')
    # Calculate bin edges using NumPy's 'auto' method
    # Calculate bin width
    bin_width = bin_edges[1] - bin_edges[0]
    set_image_para()
    sns.histplot(arr, binwidth=bin_width,color=palette1,ax = ax, kde=True)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_title(title)
    
def plot_bin2d(andata,ax,title = '',xlab = '',ylab =''):
    palette1 = sns.color_palette("colorblind",10)
    ax.scatter(andata.obs['total_counts'],andata.obs['n_genes_by_counts'], alpha=0.6,color = palette1[0],edgecolor='black')
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_title(title)
    
def plot_spatial(andata,ax,title = '',xlab = '',ylab =''):
    pass