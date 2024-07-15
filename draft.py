from utils.vizium.analysis.cellPhoneDBPrep import cellPhoneDBPrep
import os

def run(path,pathout,FilePrefix,voygerpyRead):
    # report spatial analysis vg
    SpatRep = cellPhoneDBPrep(min_counts = 1000,min_genes = 1000,min_gene_count = 1000,max_gene_count = 6281,path = path,outPath = pathout,FilePrefix = FilePrefix, voygerpyRead = voygerpyRead)
    SpatRep.save_metadata(os.path.join(pathout,'test_meta.pkl'))
    SpatRep.save_counts(os.path.join(pathout,'test_counts.pkl'))
    # SpatRep.vg_pca_report()
    # SpatRep.de_vg()
    # SpatRep.plot_vp_qc()
    # SpatRep.moransIplot_vp()
    # SpatRep.report_rank_genes_groups()
    # SpatRep.plot_nhood_enrichment()
    # SpatRep.plot_spatial_clasterList()
    # SpatRep.plot_co_occurrence()
    # SpatRep.plot_top_autocorr()
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
    FilePrefix = "_071124" 
    run(path = path, pathout = pathout, FilePrefix = FilePrefix, voygerpyRead = True)
    