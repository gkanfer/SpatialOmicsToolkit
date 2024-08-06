import rapids_singlecell as rsc
import numpy as np
import pandas as pd
import scipy.stats as stats
import cupy
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests
from matplotlib.ticker import MaxNLocator
import seaborn as sns

def rank_genes_groups(andata,return_andata = False):
    andata.X = andata.layers["log"]
    rsc.tl.rank_genes_groups_logreg(andata, groupby="cluster",groups='all')
    def z_to_p(z):
        return 2 * (1 - stats.norm.cdf(abs(z)))
    p_values = pd.DataFrame.from_records(andata.uns['rank_genes_groups']['scores']).applymap(z_to_p)
    p_values_flat = p_values.values.flatten()
    _, pvals_corrected, _, _ = multipletests(p_values_flat, alpha=0.05, method='fdr_bh')
    pvals_corrected_df = pd.DataFrame(pvals_corrected.reshape(p_values.shape), columns=p_values.columns, index=p_values.index)
    if return_andata:
        andata.uns['rank_genes_groups']['p_values'] = p_values.to_numpy()
        andata.uns['rank_genes_groups']['pvals_corrected'] = pvals_corrected_df.to_numpy()
        return andata , _
    else:
        return p_values,pvals_corrected_df
    
def return_markers(andata,marker):
    if not 'rank_genes_groups' in [key for key in andata.uns.keys()]:
        andata.X = andata.layers["log"]
        rsc.tl.rank_genes_groups_logreg(andata, groupby="log",groups='all')
    def sparseTonp(data_matrix):
        if isinstance(data_matrix, (np.ndarray, np.matrix)):
            dense_data = np.array(data_matrix)
        else:
            dense_data = data_matrix.toarray()
        return dense_data
    arr = sparseTonp(andata[:,andata.var.index.isin([marker])].X)
    arr = np.ravel(arr)
    arr = cupy.float32(arr.get())
    return arr

def find_markers(andata,prob = "pvals_corrected", threshold = 0.1):
    """Identify markers based on certain criteria for each cluster in the input data.

        Args:
            andata: Anndata object containing data for marker identification.
            prob (str): The key to access p-values data in the Anndata object (default is "pvals_corrected" and "p_values").
            threshold (float): The threshold value for p-values to consider a marker (default is 0.1).

        Returns:
            dict: A dictionary where each key represents a cluster and the corresponding value is a list of markers.
    """
    scores_df = pd.DataFrame.from_records(andata.uns['rank_genes_groups']['scores'])
    names_df = pd.DataFrame.from_records(andata.uns['rank_genes_groups']['names'])
    pvalues_df = pd.DataFrame.from_records(andata.uns['rank_genes_groups'][prob])
    def return_column(vec_names,vec_score,vec_pv,hits_container,cluster,pv = pv):
        df = pd.DataFrame({"names":vec_names,"score":vec_score,"pv":vec_pv})
        arr_names_sel = df.loc[(df["score"]>0) & (df["pv"]< pv),"names"].to_list()
        if len(arr_names_sel) > 0:
            if isinstance(arr_names_sel,list):
                hits_container[cluster] = arr_names_sel
            else:
                hits_container[cluster] = arr_names_sel
            return hits_container
        else:
            return "no hits were found"
    hits_container = {}
    for cluster in range(len(names_df.columns)):
        container = return_column(vec_names = names_df.iloc[:,cluster].values, vec_score = scores_df.iloc[:,cluster].values, vec_pv = pvalues_df.iloc[:,cluster].values, hits_container = hits_container, cluster = cluster, pv = threshold)
    return container

def compute_spatial_lag(andata):
    pass
    
    
    
    