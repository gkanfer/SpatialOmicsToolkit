from utils.vizium.analysis.SpataStatReport import SpataStatReport
# from utils.vizium.analysis import analysis_sq

# def run(path,pathout,FilePrefix):
#     x = analysis_sq(path = path,outPath = pathout,FilePrefix = FilePrefix)

def run(path,pathout,FilePrefix):
    # report spatial statics
    SpatRep = SpataStatReport(path = path,outPath = pathout,FilePrefix = FilePrefix)
    # SpatRep.report_rank_genes_groups()
    SpatRep.plot_spatial()
    SpatRep.plot_nhood_enrichment()
    SpatRep.plot_spatial_clasterList()
    SpatRep.plot_co_occurrence()
    SpatRep.plot_top_autocorr()
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
    FilePrefix = "_061724" 
    run(path = path,pathout = pathout,FilePrefix = FilePrefix)
    
    
    