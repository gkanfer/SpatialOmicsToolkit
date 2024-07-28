import anndata
import pandas as pd
import pickle
from spatialomicstoolkit.report.SpataStatReport import SpataStatReport

class cellPhoneDBPrep(SpataStatReport):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def save_metadata(self, filename):
        # Extract metadata
        metadata = self.andata.obs[['cluster']].copy()
        metadata['Cell'] = metadata.index
        metadata.rename(columns={'cluster': 'cell_type'}, inplace=True)
        metadata = metadata[['Cell', 'cell_type']]
        
        # Save as pickle
        with open(filename, 'wb') as f:
            pickle.dump(metadata, f)
        print(f"Metadata file saved as {filename}")

    def save_counts(self, filename):
        # Transpose the matrix so rows are genes and columns are cells
        counts = pd.DataFrame(self.andata.X.T.todense(), index=self.andata.var_names, columns=self.andata.obs_names)
        
        # Save as pickle
        with open(filename, 'wb') as f:
            pickle.dump(counts, f)
        print(f"Counts file saved as {filename}")

# Example usage:
# adata = anndata.read_h5ad('path_to_your_anndata_file.h5ad')
# cpdb_prep = CellPhoneDBPrep(adata)
# cpdb_prep.save_metadata('test_meta.pkl')
# cpdb_prep.save_counts('test_counts.pkl')
# cellphonedb method statistical_analysis test_meta.pkl test_counts.pkl
# cellphonedb method analysis test_meta.pkl test_counts.pkl
# cellphonedb plot dot_plot --rows in/rows.txt --columns in/columns.txt
# cellphonedb plot heatmap_plot meta_data