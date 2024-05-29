from utils.vizium.addQCplots import addQCplots

def run(path, pathout, qcFilePrefixAddQCreport = 'test', qc_add=["qc_metrics"]):
    additional_plots = addQCplots(qc_add=qc_add, path=path, outPath=pathout, qcFilePrefixAddQCreport=qcFilePrefixAddQCreport)
    additional_plots.add_additional_qc_plots()

if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out"
    run(path, pathout)
    #run(path, pathout, qc_add=["qc_metrics", "rm_high_variable_genes", "quality_control_pca"])

