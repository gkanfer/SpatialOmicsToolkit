import debugpy
import scanpy as sc
import squidpy as sq
import pandas as pd
from collections import OrderedDict
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
    def __init__(self, path, outPath, FilePrefix):
        self.path = path
        self.outPath = outPath
        self.FilePrefix = FilePrefix
        self.parquet_to_csv()
        self.andata = self.readVizHD()
        
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
    
    def set_image_para(self):
        plt.rcParams['figure.dpi'] = 150
        plt.rcParams['font.family'] = ['serif']
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['axes.titlesize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        
    
    @staticmethod
    def print_all_keys(d, parent_key=''):
        """
        Prints all keys of a dictionary recursively.

        Args:
            d (dict): The dictionary to print keys from.
            parent_key (str): The parent key for nested dictionaries.

        Returns:
            None
        
        Examples:
            print_all_keys(self.andata.uns)
        """
        if isinstance(d, (dict, OrderedDict)):
            for key, value in d.items():
                new_key = f"{parent_key}.{key}" if parent_key else key
                print(new_key)
                if isinstance(value, (dict, OrderedDict)):
                    print_all_keys(value, new_key)
        
    
    
        
        
             