import rapids_singlecell as rsc
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
def rank_genes_groups(andata,return_andata = False):
    andata.X = andata.layers["counts"]
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
    