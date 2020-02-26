mcPlotMarkers <- function (seurat_obj, marker_table, n_genes = 100) {
  if (DefaultAssay(seurat_obj) != "RNA") {
    stop("Default assay is not RNA")
  }
  
  dir.create(figurePath("top-cluster-markers/"))
  
  cell_names <- as.character(unique(Idents(seurat_obj)))
  n_clust <- seq_along(unique(Idents(seurat_obj)))
  
  print("active identities:")
  print(cell_names)

  gene_table <- table(marker_table$cluster)
  seq_nums <- seq(1, n_genes, by = 20)

  parallel::mclapply(n_clust, mc.cores = 10,
    function(i) {
      genes <- marker_table[marker_table$cluster == cell_names[i],
        "Gene.name.uniq"]
      
      for(j in 1:length(seq_nums)) {
          to_plot <- genes[seq_nums[j]:(seq_nums[j] + 19)]
          if (NA %in% to_plot) {break}

          f <- FeaturePlot(seurat_obj, to_plot,
            reduction = "umap", pt.size = 0.25, combine = FALSE)
          for(k in 1:length(f)) {
            f[[k]] <- f[[k]] + NoLegend() + NoAxes()
          }

          path <- figurePath(paste0(
            "top-cluster-markers/", cell_names[i], "_top_",
            seq_nums[j], "-", (seq_nums[j] + 19), "_features.png"))
          
          png(path, width = 30, height = 25, units = "in", res = 200)
          print(cowplot::plot_grid(plotlist = f))
          dev.off()
      }
      return(cell_names[i])
      print("Done")
    }
  )
}