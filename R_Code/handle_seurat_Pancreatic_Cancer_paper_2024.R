library(RcppAnnoy)
library(Seurat)
setwd("/gpfs/gsfs10/users/kanferg/Sptial_Omics/playGround/Data/Pancreatic_Cancer_paper_2024")
seurat_st<-readRDS('PDAC_Updated.rds')
# head(seurat_st@meta.data.colnames())
# seurat_st@meta.data$Sample_ID2
# names(seurat_st@images)
seurat_list <- SplitObject(seurat_st, split.by = "Sample_ID2")

##### Test 
sample_name <- names(seurat_list)[1]
seurat_object <- seurat_list[[sample_name]]
Assays(seurat_object)
origen <- as.data.frame(seurat_object@meta.data$Origin)
origen_text <- unique(origen)
paste(origen_text, paste0(sample_name, "_raw_counts.csv"),sep = "#")
raw_counts <- as.data.frame(GetAssayData(seurat_object, assay = "Spatial", slot = "counts"))
xy_coords <- GetTissueCoordinates(seurat_object@images$IU_PDA_HM9)
as.data.frame(seurat_object@meta.data$Origin)
######

# saving data to csv

for (image_name in names(seurat_list)){
  print(image_name)
  seurat_subset <- seurat_list[[image_name]]
  raw_counts <- GetAssayData(seurat_subset, assay = "Spatial", slot = "counts")
  image_id = names(seurat_subset@images)
  xy_coords <- as.data.frame(GetTissueCoordinates(seurat_subset@images[[image_id]]))
  origen <- as.data.frame(seurat_subset@meta.data$Origin)
  origen_text <- unique(origen)
  cell_type <- as.data.frame(seurat_subset$first_type)
  clusters <- as.data.frame(seurat_subset$CompositionCluster_CC)
  var <- rownames(raw_counts)
  # Save the data to CSV files
  write.csv(var,paste0(image_id, "_var.csv"))
  write.csv(as.matrix(raw_counts),paste0(image_id, "_raw_counts.csv"))
  write.csv(xy_coords,paste0(image_id, "_xy_coordinates.csv"))
  write.csv(origen, paste0(image_id, "_origen.csv"))
  write.csv(cell_type, paste0(image_id, "_cell_type.csv"))
  write.csv(clusters, paste0(image_id, "_clusters.csv"))
  
  # write.csv(as.matrix(raw_counts),paste0(image_name, "_raw_counts.csv"))
  # write.csv(xy_coords,paste0(image_name, "_xy_coordinates.csv"))
  # write.csv(origen, paste0(image_name, "_origen.csv"))
  # write.csv(cell_type, paste0(image_name, "_cell_type.csv"))
  # write.csv(clusters, paste0(image_name, "_clusters.csv"))
  
  # write.csv(as.matrix(raw_counts),paste(origen_text,paste0(image_name, "_raw_counts.csv"),sep="#"))
  # write.csv(xy_coords,paste(origen_text,paste0(image_name, "_xy_coordinates.csv"),sep="#"))
  # write.csv(origen, paste(origen_text,paste0(image_name, "_origen.csv"),sep="#"))
  # write.csv(cell_type, paste(origen_text,paste0(image_name, "_cell_type.csv"),sep="#"))
  # write.csv(clusters, paste(origen_text,paste0(image_name, "_clusters.csv"),sep="#"))
}





