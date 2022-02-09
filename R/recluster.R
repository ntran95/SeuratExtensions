# recluster dimensions reduction on subsets of Seurat object 
recluster_subsets <- function(seurat_obj,dims,res,path){
  seurat_obj <- NormalizeData(
    seurat_obj, verbose = FALSE)
  
  seurat_obj <- FindVariableFeatures(
    seurat_obj, selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
  
  seurat_obj <- ScaleData(seurat_obj, 
                                    verbose = TRUE)
  
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, 
                                 verbose = TRUE)
  
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, 
                                        k.param = 20)
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  
  seurat_obj <- BuildClusterTree(
    seurat_obj, reorder = TRUE, reorder.numeric = TRUE )
  
  saveRDS(object = seurat_obj, file = path)
  
  return(seurat_obj)
  
}

