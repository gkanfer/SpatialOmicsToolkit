from utils.viziumHD import viziumHD

def run(path,pathout,plotPDF = ["qc","embd"]):
    if plotPDF == "qc":
        viziumHD(path = path,outPath = pathout,totalThr = 100,geneThr = 100 , bins_gene_plot = 20,bins_gene_plot_thr = 20,qcFilePrefix = '008um_zoom_thr_100')
    elif plotPDF == "embd":
        viziumHD(path = path,outPath = pathout,totalThr = 1500,geneThr = 1500 , bins_gene_plot = 20,bins_gene_plot_thr = 20,qcFilePrefix = '008um',qc = False).filterANDnorm(max_counts = 1500, max_genes = 100)
    
if __name__ == "__main__":
    path = "/data/kanferg/Sptial_Omics/playGround/Data/Visium_HD_Mouse_Brain_square_example/square_008um"
    pathout = "/data/kanferg/Sptial_Omics/playGround/out"
    run(path,pathout,plotPDF = 'qc')
    
    
    

'''
self.totalThr = totalThr
self.bins_total = bins_total
self.bins_gene_plot = 200
self.geneThr = geneThr∏
self.bins_gene_plot_thr = bins_gene_plot_thr
self.qcFilePrefix = qcFilePrefix
'''