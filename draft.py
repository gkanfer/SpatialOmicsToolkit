def agrrData(file_select):
    def get_file(attrib,file_sel):
        file_ind = [file for file in file_sel if attrib in file]
        return file_ind[0]
    raw_counts_ind =  get_file(attrib = 'raw_counts',file_sel = file_select)   
    xy_coordinates_ind = get_file(attrib = 'xy_coordinates',file_sel =file_select)
    var_ind = get_file(attrib = 'var',file_sel =file_select)
    cell_type_ind = get_file(attrib = 'cell_type',file_sel =file_select)
    clusters_ind = get_file(attrib = 'clusters',file_sel =file_select)
    # read data
    counts_df = pd.read_csv(os.path.join(path_016,raw_counts_ind),index_col=None)
    counts = counts_df.iloc[:,1:].to_numpy()
    del counts_df
    xy_coordinates = pd.read_csv(os.path.join(path_016,xy_coordinates_ind),index_col=None).iloc[:,1:].to_numpy()
    var = pd.read_csv(os.path.join(path_016,var_ind),index_col=None).iloc[:,1:].rename(columns = {"x":'id'})
    obs = pd.read_csv(os.path.join(path_016,cell_type_ind),index_col=None).iloc[:,1:].rename(columns = {'seurat_subset$first_type':'cell_type'})
    obs['clusters'] = pd.read_csv(os.path.join(path_016,clusters_ind),index_col=None).iloc[:,1:].rename(columns = {'seurat_subset$CompositionCluster_CC':'clusters'})
    andata = AnnData(counts.T,var = var ,obsm={"spatial": xy_coordinates}, obs = obs, uns = {"sample_name":file_select})
    andata.obsm['spatial'] = np.array(andata.obsm['spatial'], dtype=np.float64)
    andata.var = andata.var.set_index('id')
    andata.var.index.name = None
    andata.X = sp.csr_matrix(andata.X)
    return andata