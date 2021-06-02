
testPCs <- function(seurat_obj, 
                    from_to = 5:10,
                    include_conditions = FALSE,
                    group.by = "data.set",
                    n_cores = 10, 
                    res = 1.2,
                    genes = NULL) {
  dir.create(figurePath("exploratory-PCs/"))
  dir.create(figurePath(paste0("exploratory-PCs/dims_", 
                               min(from_to), "_to_", max(from_to),
                               "_res_", res, "/")))
  PCs <- from_to
  
  parallel::mclapply(PCs, function(i){
    dims <- c(1:i)
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims, k.param = 20)
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    seurat_obj <- RunUMAP(seurat_obj, dims = dims)
    
    seurat_obj <- BuildClusterTree(seurat_obj,
                                   reorder = TRUE, reorder.numeric = TRUE)
    
    print("printing to png...")
    umap_clusters <- DimPlot(seurat_obj, reduction = "umap",
                                     pt.size = 0.20,
                             label = TRUE, label.size = 4)
    
    umap_clusters <- cleanUMAP(umap_clusters)
    
    png(figurePath(paste0("exploratory-PCs/dims_", 
                          min(from_to), "_to_", max(from_to),
                          "_res_", res, "/", "UMAP_by_clusters_",
                          "dims_", max(dims),
                          "_res_", res,".png")), 
        width = 12, height = 10, units = "in", res = 300)
    print(umap_clusters)
    dev.off()
    
    if(!is.null(genes)){
      print("plotting feature plots...")
      DefaultAssay(seurat_obj) <- 'RNA'
      e <- FeaturePlot(seurat_obj, genes,
                       reduction = "umap", pt.size = 0.25, 
                       combine = FALSE,label = TRUE,label.size = 2)
      
      for (i in 1:length(e)) {
        e[[i]] <- e[[i]] + NoLegend() + NoAxes() +
          theme(
            plot.subtitle = element_text(hjust = 0.5)
          )
      }
      
      plot_list <- cowplot::plot_grid(plotlist = e)
      
      #save_plot is a wrapper of ggsave but automatically scales plot_grid 
      #depending on input of ncol or nrow
      cowplot::save_plot(
        filename = figurePath(paste0("exploratory-PCs/dims_", 
                                     min(from_to), "_to_", max(from_to),
                                     "_res_", res, "/", "FeaturePlots", "_",
                                     "dims_", max(dims), "_res_", res,".png")),
        plot = plot_list,ncol = 4,base_asp = .25)
      
    }
    
    if(include_conditions) {
      umap_dataset <- DimPlot(seurat_obj, reduction = "umap",
                                      pt.size = 0.25,
                              label = TRUE, label.size = 4, 
                              group.by = group.by)
      
      umap_dataset <- cleanUMAP(umap_dataset)
      
      png(figurePath(paste0("exploratory-PCs/dims_", 
                            min(from_to), "_to_", max(from_to),
                            "_res_", res, "/", "UMAP_by_", group.by, "_",
                            "dims_", max(dims), "_res_", res,".png")),
          width = 14, height = 12,
          units = "in", res = 300)
      print(umap_dataset)
      dev.off()
    }
    
   # FeaturePlot of nFeature_RNA, heuristic approach to determining whether 
    #clustering is influenced by cell quality. 
    #use to determine overclustering bc communitites will aggregate by high and 
    #low nFeature_RNA
    
  f <- FeaturePlot(seurat_obj,features = "nFeature_RNA", 
                   pt.size = 0.25,
                   label = TRUE,label.size = 2) + NoAxes()
  
  png(figurePath(paste0("exploratory-PCs/dims_", 
                        min(from_to), "_to_", max(from_to),
                        "_res_", res, "/", "FeaturePlot_by_nFeature", "_",
                        "dims_", max(dims), "_res_", res,".png")),
      width = 14, height = 12,
      units = "in", res = 300)
  print(f)
  dev.off()
  
  if(save){
    saveRDS(seurat_obj, file = dataPath(paste0(
      "SeurObj_", script_name, "_dim_", max(dims), "_res_", res,"_.RDS")))
              }
  
  }, mc.cores = n_cores)
}

# Range of PCs to test
# Execute

# testPCs2 <- function(seurat_obj, from_to = 5:10,
#   include_conditions = FALSE, n_cores = 10, res = 1.2, 
#   save = FALSE) {
#   
#     meta_common_features <- read.table(
#     file = "/n/projects/nt2473/Analysis/Data/gene-lists/meta_common_features.tsv", 
#     sep = "", header = T)
#     
#     common_features <- scan(paste0("/n/projects/nt2473/Analysis/",
#                                    "Data/gene-lists/common_features.txt"),
#                             what = "character")
#     
#     dir.create(figurePath("exploratory-PCs/"),showWarnings = FALSE)
#     PCs <- from_to
#     
#     parallel::mclapply(PCs, function(i){
#       dims <- c(1:i)
#       seurat_obj <- FindNeighbors(seurat_obj, dims = dims, k.param = 20)
#       seurat_obj <- FindClusters(seurat_obj, resolution = res)
#       seurat_obj <- RunUMAP(seurat_obj, dims = dims,
#                             min.dist = 0.3, spread = .5)
# 
#       seurat_obj <- BuildClusterTree(seurat_obj,
#         reorder = TRUE, reorder.numeric = TRUE)
# 
#       umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.50,
#         label = TRUE, label.size = 4)
# 
#       umap_clusters <- cleanUMAP(umap_clusters)
# 
#       png(figurePath(paste0("exploratory-PCs/UMAP_Dims_",
#         length(dims),".png")), width = 12, height = 10, units = "in", res = 200)
#       print(umap_clusters)
#       dev.off()
#       
#       DefaultAssay(seurat_obj) <- 'RNA'
#       e <- FeaturePlot(seurat_obj, common_features,
#                        reduction = "umap", pt.size = 0.25, combine = FALSE, 
#                        label = TRUE)  
#       
#       for (i in 1:length(e)) {
#         e[[i]] <- e[[i]] + NoLegend() + NoAxes() + 
#           labs(subtitle = meta_common_features$marker.ident[i]) +
#           theme(
#             plot.subtitle = element_text(hjust = 0.5)
#           )
#       }
#       
#       png(figurePath(paste0("exploratory-PCs/common_features_Dims_", 
#         length(dims),
#       ".png")), width = 40,
#           height = 80, units = "in", res = 200)
#       print(cowplot::plot_grid(plotlist = e, ncol = 4))
#       dev.off()
# 
#       if(include_conditions) {
#         umap_dataset <- DimPlot(seurat_obj, reduction = "umap",pt.size = 0.25,
#           label = FALSE, label.size = 4, group.by = "data.set",
#           cols = trt_colors)
# 
#         umap_dataset <- cleanUMAP(umap_dataset)
# 
#         png(figurePath(paste0("exploratory-PCs/UMAP_", "Dataset",
#           length(dims), "_res", res,".png")), width = 14, height = 12,
#           units = "in", res = 200)
#         print(umap_dataset)
#         dev.off()
#       }
# 
#       if(save){
#         saveRDS(seurat_obj, file = dataPath(paste0(
#           "SeurObj_", script_name, "_dim_", max(dims),"_.RDS")))
#           }
# 
#   }, mc.cores = n_cores)
# }

 # Range of PCs to test
 # Execute
