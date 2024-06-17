from utils.vizium.qc.applyqc import applyqc

def run(path,pathout,FilePrefix):
    # report qc
    qc = applyqc(path = path,outPath = pathout,FilePrefix = FilePrefix)
    qc.qc_report()
    qc.plot_mitochondrial_data()
    qc.perform_pca_and_save_plot()
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_016um"
    pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_1"
    FilePrefix = "_061724" 
    run(path = path,pathout = pathout,FilePrefix = FilePrefix)