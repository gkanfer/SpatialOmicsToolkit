import scanpy as sc
import cupy as cp
import os
import time
import rapids_singlecell as rsc
import numpy as np
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)
import zarr
from collections import OrderedDict
from scipy.sparse import csr_matrix
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from scipy.sparse import csr_matrix
import scipy
import anndata
from collections import OrderedDict

class StDatareader_rsc:
    '''
    A class to handle the processing and quality control reporting of Visium high-definition spatial transcriptomics data.
    
    Parameters
    ----------
    path : str
        The path to the directory containing the spatial transcriptomics data files.
    outPath : str
        The output directory path where the results, including the quality control report, will be saved.
    '''
    def __init__(self, path, outPath, FilePrefix, hdffileName = "sp_countAndata.h5ad", method = "vizium", subsample = None):
        self.path = path
        self.outPath = outPath
        self.FilePrefix = FilePrefix
        self.hdffileName = hdffileName
        self.method = method
        self.subsample = subsample
        if self.method == "vizium":
            self.parquet_to_csv()
            self.andata = self.readVizHD()
        else:
            self.andata = self.readXenium()
        
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
    
    def readXenium(self):
        path_xenium = os.path.join(self.path,"cell_feature_matrix.h5")
        path_cells = os.path.join(self.path,"cells.zarr.zip")
        adata = sc.read_10x_h5(path_xenium)
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="mt-")
        rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT"])
        def open_zarr(path: str) -> zarr.Group:
            store = (zarr.ZipStore(path, mode="r") if path.endswith(".zip") else zarr.DirectoryStore(path))
            return zarr.group(store=store)
        root = open_zarr(path_cells)
        column_names = dict(root['cell_summary'].attrs.items())['column_names']
        def build_obs(andata,root,column_names):
            for i in range(len(column_names)):
                andata.obs[str(column_names[i])] = np.array(root["cell_summary"])[:,i]
            spatial = andata.obs[["cell_centroid_x", "cell_centroid_y"]]
            adata.obsm["spatial"] = spatial.values
            return andata
        return build_obs(adata,root,column_names)
    
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
    
    def writeAnn(self):
        print("prepeare for writing")
        self.andata.write(os.path.join(self.outPath, self.hdffileName))    
    
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
        
    def writeAnn_xenium(self):
        print("prepare for writing")
        self.andata.X = self.andata.layers['counts']
        del self.andata.layers
        keys =[keys for keys in self.andata.obsp.keys()]
        for key in keys:
            matrix = self.andata.obsp[key]
            if isinstance(matrix, scipy.sparse.spmatrix):
                self.andata.obsp[key] = np.array(self.andata.obsp[key].todense())
            if isinstance(matrix, OrderedDict):
                self.andata.obsp[key] = dict(self.andata.obsp[key]) 
                # self.andata.obsp.pop(key)
        keys = [keys for keys in self.andata.obsm.keys()]
        for key in keys:
            matrix = self.andata.obsm[key]
            if isinstance(matrix, scipy.sparse.spmatrix):
                self.andata.obsm[key] = np.array(self.andata.obsm[key].todense())
            if isinstance(matrix, OrderedDict):
                self.andata.obsm[key] = dict(self.andata.obsm[key]) 
                # self.andata.obsm.pop(key)
        keys = [keys for keys in self.andata.uns.keys()]
        for key in keys:
            matrix = self.andata.uns[key]
            if isinstance(matrix, scipy.sparse.spmatrix):
                self.andata.uns[key] = np.array(self.andata.uns[key].todense()) 
                #self.andata.uns.pop(key)
            if isinstance(matrix, OrderedDict):
                self.andata.uns[key] = dict(self.andata.uns[key])
        del self.andata.uns['config']
        del self.andata.uns['spatial']['knn_weights']
        self.andata.obsm['geometry'].to_parquet(path = os.path.join(self.outPath,self.FilePrefix + 'test.parquet'), compression="snappy")
        del self.andata.obsm['geometry'] 
        self.andata.write(os.path.join(self.outPath, self.hdffileName))
        
        
             