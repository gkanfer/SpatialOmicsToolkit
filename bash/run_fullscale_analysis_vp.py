# Explanation of Parameters and Workflow
# In the small-scale run, I used the following parameters to reduce the application run time:

# min_counts = 1000
# min_genes = 1000
# min_gene_count = 1000
# max_gene_count = 6281
# These settings filtered the data, speeding up the process. The results were saved in the out2 directory.

# For testing the new VoyagerPy functions from scratch, I applied a more relaxed filtering strategy, which resulted in a much longer running time. Therefore, I used sbatch to run the analysis efficiently.

from utils.vizium.analysis.SpataStatReport import SpataStatReport

def run(path,pathout,FilePrefix,voygerpyRead):
    # report spatial analysis vg
    SpatRep = SpataStatReport(min_counts = 100, min_genes = 80,path = path,outPath = pathout,FilePrefix = FilePrefix, voygerpyRead = voygerpyRead)
    # save file
    SpatRep.writeAnn()
    # reports 
    SpatRep.vg_pca_report()
    SpatRep.plot_spatial_vp()
    SpatRep.de_vg()
    SpatRep.plot_vp_qc()
    SpatRep.moransIplot_vp()
    SpatRep.hist_2d_plot()
    
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_2"
    FilePrefix = "_071324" 
    run(path = path, pathout = pathout, FilePrefix = FilePrefix, voygerpyRead = True)