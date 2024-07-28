from spatialomicstoolkit.analysis.cellPhoneDBPrep import cellPhoneDBPrep

def run(path,pathout,FilePrefix,voygerpyRead):
    # report spatial analysis vg
    SpatRep = cellPhoneDBPrep(path = path, outPath = pathout, FilePrefix = FilePrefix, voygerpyRead = voygerpyRead, method = "xenium", subsample = False, min_counts=10,min_genes = 10,min_gene_count = 3)
    #report
    SpatRep.vg_pca_report()
    SpatRep.de_vg()
    SpatRep.moransIplot_vp()
    SpatRep.plot_vp_qc()
    SpatRep.writeAnn_xenium()
    #cellPhone
    SpatRep.save_metadata('/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1/072223_meta.pickle')
    SpatRep.save_counts('/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1/072223_counts.pickle')
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
    FilePrefix = "_072224" 
    run(path = path, pathout = pathout, FilePrefix = FilePrefix, voygerpyRead = True)