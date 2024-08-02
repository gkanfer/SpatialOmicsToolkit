import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import seaborn as sns
import os
import gzip
import numpy as np
import rapids_singlecell as rsc
import pandas as pd
import random


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
        arr = andata.obs[column].values
    else:
        arr = andata.var[column].values
    bin_edges = np.histogram_bin_edges(arr, bins='auto')
    # Calculate bin edges using NumPy's 'auto' method
    # Calculate bin width
    bin_width = bin_edges[1] - bin_edges[0]
    set_image_para()
    sns.histplot(arr, binwidth=bin_width,color=palette1,ax = ax, kde=True)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_title(title)
    
def custom_paramsForSPatialPlot():
    custom_params = {"xtick.labelsize": 0,      
                    "ytick.labelsize": 0,      
                    "axes.labelsize": 0,       
                    "xtick.major.size": 0,     
                    "xtick.minor.size": 0,     
                    "ytick.major.size": 0,    
                    "ytick.minor.size": 0 }
    return custom_params
    
def plot_bin2d(andata,ax,title = '',xlab = '',ylab =''):
    palette1 = sns.color_palette("colorblind",10)
    ax.scatter(andata.obs['total_counts'],andata.obs['n_genes_by_counts'], alpha=0.6,color = palette1[0],edgecolor='black')
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_title(title)
    
def plot_spatial_qc(andata,ax,column,title = '',xlab = '',ylab ='',size = 2):
    '''
    e.g.
    # Create a blue gradient palette
    palette = sns.color_palette("Blues", as_cmap=True)
    listed_cmap = ListedColormap(palette(np.linspace(0, 1, 256)))

    # Plot using Matplotlib
    fig, ax = plt.subplots(figsize=(5, 5))
    sc = ax.scatter(x=df['x'], y=df['y'], c=df['total_counts'], cmap=listed_cmap)

    # Add color bar with title
    cbar = fig.colorbar(sc)
    #cbar.set_label('Total Count', fontsize=12)

    # Automatically set tick locations and labels
    locator = MaxNLocator(nbins=5)
    cbar.locator = locator
    cbar.update_ticks()

    # Explicitly set tick labels to ensure they show
    ticks = cbar.get_ticks()
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{int(tick)}' for tick in ticks])

    # Remove ticks and labels from the main plot
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.show()
    '''
    df = pd.DataFrame({str(column):andata.obs[column],'x':andata.obsm['spatial'][:,0],'y':andata.obsm['spatial'][:,1]})
    custom_params = custom_paramsForSPatialPlot()
    sns.set_theme(style="whitegrid", palette="pastel", rc=custom_params)
    palette = sns.color_palette("Blues", as_cmap=True)
    listed_cmap = ListedColormap(palette(np.linspace(0, 1, 256)))
    sc = ax.scatter(x=df['x'], y=df['y'], c=df['total_counts'], cmap=listed_cmap,s = size)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    return sc

def plot_spatial(andata,ax,features = None,title = '',xlab = '',ylab ='',size = 2, random_palette = False):
    palette = sns.color_palette("tab20") + sns.color_palette("tab20b") + sns.color_palette("tab20c")
    if random_palette:
        random.shuffle(palette) 
    df = pd.DataFrame({'cluster':andata.obs['cluster'],'x':andata.obsm['spatial'][:,0],'y':andata.obsm['spatial'][:,1]})
    if features:
        df[df['cluster'].isin([features])]
    num_classes = len(df['cluster'].unique())
    if num_classes==1:
        listed_cmap = ListedColormap(palette)
    else:
        num_classes = len(df['cluster'].unique())
        extended_palette = palette * (num_classes // len(palette) + 1)
        extended_palette = extended_palette[:num_classes]
        listed_cmap = ListedColormap(extended_palette)
        
    custom_params = custom_paramsForSPatialPlot()
    sns.set_theme(style="whitegrid", palette="pastel", rc=custom_params)
    for i, cluster in enumerate(df['cluster'].unique()):
        cluster_data = df[df['cluster'] == cluster]
        ax.scatter( x=cluster_data['x'], y=cluster_data['y'], color=listed_cmap(i), label=f'{cluster}', s=1, alpha=0.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    legend = ax.legend( title="Cluster",
                        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
                        loc='upper left',
                        fontsize='small',  # Control the font size
                        title_fontsize='medium',
                        markerscale=5,  # Increase the size of the legend markers
                        frameon=False# Control the title font size
                        )
    
def plot_expression(df,marker,ax,**kwargs):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)  
    g = sns.violinplot(data = df, x = "cluster", y = marker, ax = ax, **kwargs)
    max_value = df[marker].max()
    g.set(ylim = (0,max_value+1))
    ax.text(0.95, 0.95, marker, transform=ax.transAxes, fontsize=12, verticalalignment='top',horizontalalignment='right')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_ylabel('read count \n [log-normalized]')
    ax.set_xlabel('')