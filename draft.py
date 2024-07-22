from spatialomicstoolkit.report.SpataStatReport import SpataStatReport

def run(path,pathout,FilePrefix,voygerpyRead):
    # report spatial analysis vg
    SpatRep = SpataStatReport(path = path, outPath = pathout, FilePrefix = FilePrefix, voygerpyRead = voygerpyRead, method = "xenium", subsample = True)
    SpatRep.vg_pca_report()
    SpatRep.de_vg()
    SpatRep.moransIplot_vp()
    SpatRep.plot_vp_qc()
    # SpatRep.report_rank_genes_groups()
    # SpatRep.plot_nhood_enrichment()
    # SpatRep.plot_spatial_clasterList()
    # SpatRep.plot_co_occurrence()
    # SpatRep.plot_top_autocorr()
    
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Xenium/output_temp"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
    FilePrefix = "_072124" 
    run(path = path, pathout = pathout, FilePrefix = FilePrefix, voygerpyRead = True)