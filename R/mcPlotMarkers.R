mcPlotMarkers <- function (seurat_obj, marker_table, n_genes = 100) {
  dir.create(figurePath("top-cluster-markers/"))
  all_BCs <- Idents(seurat_obj)
  
  cell_names <- as.character(unique(Idents(seurat_obj)))
  n_clust <- seq_along(unique(Idents(seurat_obj)))
  
  print("active identities:")
  print(cell_names)

  if ("cell.type.ident" %in% colnames(marker_table))

  gene_table <- table(cluster_column)
  seq_nums <- seq(1, n_genes, by = 20)

  parallel::mclapply(n_clust, mc.cores = 10,
    function(i) {
      genes <- marker_table[cluster_column == cell_names[i],
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
            seq_nums[j], "-", (seq_nums[j] + 19), "_features2.png"))
          png(path, width = 30, height = 25, units = "in", res = 200)
          print(cowplot::plot_grid(plotlist = f))
          dev.off()
      }
    }
  )
}