
testPCs <- function(seurat_obj, from_to = 5:10,
  include_conditions = FALSE, n_cores = 10) {
    dir.create(figurePath("exploratory-PCs/"))
    PCs <- from_to
    
    parallel::mclapply(PCs, function(i){
      dims <- c(1:i)
      seurat_obj <- FindNeighbors(seurat_obj, dims = dims, k.param = 20)
      seurat_obj <- FindClusters(seurat_obj, resolution = 1.0)
      seurat_obj <- RunUMAP(seurat_obj, dims = dims)

      seurat_obj <- BuildClusterTree(seurat_obj,
        reorder = TRUE, reorder.numeric = TRUE)

      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.20,
        label = TRUE, label.size = 4)

      umap_clusters <- cleanUMAP(umap_clusters)

      png(figurePath(paste0("exploratory-PCs/UMAP_Clusters",
        length(dims),".png")), width = 12, height = 10, units = "in", res = 200)
      print(umap_clusters)
      dev.off()

      if(include_conditions) {
        umap_dataset <- DimPlot(seurat_obj, reduction = "umap",pt.size = 0.25,
          label = FALSE, label.size = 4, group.by = "data.set",
          cols = trt_colors)

        png(figurePath(paste0("exploratory-PCs/UMAP_", "Dataset",
          length(dims),".png")), width = 14, height = 12,
          units = "in", res = 200)
        print(umap_dataset)
        dev.off()
      }
    }, mc.cores = n_cores)
  }

 # Range of PCs to test
 # Execute


