from spatialomicstoolkit.analysis.SpatialStatsSQ import SpatialStatsSQ
import os

def run(path,pathout,FilePrefix,voygerpyRead):
    # report spatial analysis vg
    SpatRep = SpatialStatsSQ(path = path, outPath = pathout, FilePrefix = FilePrefix, voygerpyRead = voygerpyRead, method = "xenium", subsample = 100_000, min_counts = 10, min_genes = 10, min_gene_count = 3)
    SpatRep.perform_analysis_vp()
    # SpatRep.vg_pca_report()
    # SpatRep.vg_pca_report()
    # SpatRep.de_vg()
    # SpatRep.moransIplot_vp()
    # SpatRep.plot_vp_qc()
    SpatRep.writeAnn_xenium()
    # SpatRep.save_metadata('/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2/test_meta.pkl')
    # SpatRep.save_counts('/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2/test_counts.pkl')
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
    FilePrefix = "_072324" 
    run(path = path, pathout = pathout, FilePrefix = FilePrefix, voygerpyRead = True)