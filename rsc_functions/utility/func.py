


def writeAnn_xenium(andata,outPath,hdffileName,FilePrefix):
    print("prepare for writing")
    andata.X = andata.layers['counts']
    del andata.layers
    keys = [keys for keys in andata.obsp.keys()]
    for key in keys:
        matrix = andata.obsp[key]
        if isinstance(matrix, scipy.sparse.spmatrix):
            andata.obsp[key] = np.array(andata.obsp[key].todense())
        if isinstance(matrix, OrderedDict):
            andata.obsp[key] = dict(andata.obsp[key]) 
            # andata.obsp.pop(key)
    keys = [keys for keys in andata.obsm.keys()]
    for key in keys:
        matrix = andata.obsm[key]
        if isinstance(matrix, scipy.sparse.spmatrix):
            andata.obsm[key] = np.array(andata.obsm[key].todense())
        if isinstance(matrix, OrderedDict):
            andata.obsm[key] = dict(andata.obsm[key]) 
            # andata.obsm.pop(key)
    keys = [keys for keys in andata.uns.keys()]
    for key in keys:
        matrix = andata.uns[key]
        if isinstance(matrix, scipy.sparse.spmatrix):
            andata.uns[key] = np.array(andata.uns[key].todense()) 
            #andata.uns.pop(key)
        if isinstance(matrix, OrderedDict):
            andata.uns[key] = dict(andata.uns[key])
    del andata.uns['config']
    del andata.uns['spatial']['knn_weights']
    andata.obsm['geometry'].to_parquet(path = os.path.join(outPath,FilePrefix + 'test.parquet'), compression="snappy")
    del andata.obsm['geometry'] 
    andata.write(os.path.join(outPath, hdffileName))