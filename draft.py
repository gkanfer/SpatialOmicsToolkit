import STAligner
import os

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg

import scipy
import networkx

import torch

used_device = torch.device('cuda')

pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_4"
path_age_58 = "/data/kanferg/Sptial_Omics/playGround/Data/Breast_Cancer/age_58/binned_outputs/square_016um"
path_age_76 = "/data/kanferg/Sptial_Omics/playGround/Data/Breast_Cancer/age_76/binned_outputs/square_016um"
path_list = [path_age_58,path_age_76]

def read_data(path_data):
    def parquet_to_csv(path):
        '''
        Converts a Parquet file to a CSV file if the CSV file does not already exist.
        '''
        file_path = os.path.join(path,'spatial/tissue_positions_list.csv')
        if not os.path.exists(file_path):
            df = pd.read_parquet(os.path.join(path,'spatial/tissue_positions.parquet'))
            # Write to a CSV file
            df.to_csv(os.path.join(path,'spatial/tissue_positions_list.csv'), index=False)
        return
    parquet_to_csv(path_data)
    andata = sc.read_visium(path=path_data,load_images=False)
    positions = pd.read_csv(os.path.join(path_data,'spatial/tissue_positions_list.csv'),index_col=0,)
    positions.columns = [
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_col_in_fullres",
                "pxl_row_in_fullres",
            ]
    andata.obs = andata.obs.join(positions, how="left")
    andata.obsm["spatial"] = andata.obs[
                ["pxl_row_in_fullres", "pxl_col_in_fullres"]
            ].to_numpy()
    andata.obs.drop(
        columns=["pxl_row_in_fullres", "pxl_col_in_fullres"],
        inplace=True,
    )
    andata.obsm['spatial'] = np.array(andata.obsm['spatial'], dtype=np.float64)
    andata.var_names_make_unique()
    return andata

Batch_list = []
adj_list = []
section_ids = ['0','1']
for i,p in enumerate(path_list):
    adata = read_data(p)
    adata.obs_names = [x + '_' + section_ids[i] for x in adata.obs_names]
    STAligner.Cal_Spatial_Net(adata, rad_cutoff=1.3)  
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,flavor = 'cell_ranger',n_top_genes = 300)
    
    adata = adata[:, adata.var['highly_variable']]
    adj_list.append(adata.uns['adj'])
    Batch_list.append(adata)




adata_concat = ad.concat(Batch_list, label="slice_name", keys=section_ids)
adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')


##### Construct edgeList using memory efficient function
# adj_concat = np.asarray(adj_list[0].todense())
# adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[1].todense()))
# adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[2].todense()))
# adata.uns['edgeList'] = np.nonzero(adj_concat)

def memory_efficient_block_diag(adj_list):
    """
    Constructs a block diagonal matrix from a list of sparse adjacency matrices in a memory-efficient way.
    
    Parameters:
        adj_list (list): A list of sparse matrices (scipy.sparse).

    Returns:
        scipy.sparse.csr_matrix: Block diagonal sparse matrix.
    """
    # Calculate the total size of the final matrix
    total_rows = sum(adj.shape[0] for adj in adj_list)
    total_cols = sum(adj.shape[1] for adj in adj_list)

    # Initialize sparse matrices for the data, row indices, and column indices
    data = []
    row_indices = []
    col_indices = []

    row_offset = 0
    col_offset = 0

    for adj in adj_list:
        # Convert sparse matrix to COO format for efficient indexing
        adj_coo = adj.tocoo()

        # Append the data and adjusted indices
        data.append(adj_coo.data)
        row_indices.append(adj_coo.row + row_offset)
        col_indices.append(adj_coo.col + col_offset)

        # Update offsets for the next block
        row_offset += adj.shape[0]
        col_offset += adj.shape[1]

    # Concatenate data and indices from all blocks
    data = np.concatenate(data)
    row_indices = np.concatenate(row_indices)
    col_indices = np.concatenate(col_indices)

    # Construct the final sparse matrix
    block_diag_matrix = sp.csr_matrix((data, (row_indices, col_indices)), shape=(total_rows, total_cols))
    return block_diag_matrix
# Example usage
adj_concat = memory_efficient_block_diag(adj_list)
adata_concat.uns['edgeList'] = np.nonzero(adj_concat)

iter_comb = [(1, 0)] ## Fix slice 0 as reference to align
adata_concat = STAligner.train_STAligner(adata_concat, verbose=True, knn_neigh = 10, iter_comb = iter_comb,
                                                        margin=2.5,  device=used_device)

