# from local mac

# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
library(Seurat)
library(SeuratData)
library(SeuratDisk)
# devtools::install_github("satijalab/seurat-data")
# devtools::install_github("mojaveazure/seurat-disk")
setwd('/Users/kanferg/Desktop/Projects/Sptial_Omics_playGround/data/Pancreatic_Cancer_paper_2024')
seurat_st<-readRDS('/Users/kanferg/Desktop/Projects/Sptial_Omics_playGround/data/Pancreatic_Cancer_paper_2024/PDAC_Updated.rds')
SaveH5Seurat(seurat_st, filename = "PDAC_Updated.h5ad")
Convert("PDAC_Updated.h5ad", dest = "h5ad")