import debugpy
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import seaborn as sns
import os
import gzip
import numpy as np

class viziumHD:
    '''
    A class to handle the processing and quality control reporting of Visium high-definition spatial transcriptomics data.
    
    Parameters
    ----------
    path : str
        The path to the directory containing the spatial transcriptomics data files.
    outPath : str
        The output directory path where the results, including the quality control report, will be saved.
    '''
    def __init__(self, path, outPath):
        self.path = path
        self.outPath = outPath
        self.parquet_to_csv()
        self.andata = self.readVizHD()
        plt.rcParams['figure.dpi'] = 150
        plt.rcParams['font.family'] = ['serif']
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['axes.titlesize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        
        
    def parquet_to_csv(self):
        '''
        Converts a Parquet file to a CSV file if the CSV file does not already exist.
        '''
        file_path = os.path.join(self.path,'spatial/tissue_positions_list.csv')
        if not os.path.exists(file_path):
            df = pd.read_parquet(os.path.join(self.path,'spatial/tissue_positions.parquet'))
            # Write to a CSV file
            df.to_csv(os.path.join(self.path,'spatial/tissue_positions_list.csv'), index=False)
        return
        
    def readVizHD(self):
        return sc.read_visium(path = self.path)
    
    def printAnnD(self):
        print(f'{self.andata}')
        return self
    
    
        
        
             