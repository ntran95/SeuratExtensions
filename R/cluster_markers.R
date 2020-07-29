mcFindMarkers <- function(seurat_obj, file_prefix = "", save_table = TRUE,
  save_raw = TRUE, pval_cutoff = 0.05, n_cores = 10) {
  
  if (DefaultAssay(seurat_obj) != "RNA") {
  stop("Default assay is not RNA")
  }

  cell_names <- unique(Idents(seurat_obj))
  n_clust <- seq_along(unique(Idents(seurat_obj)))

  print("Clusters to evaluate:")
  print(cell_names)

  marker_results <- list()[n_clust]
  marker_results <- parallel::mclapply(
    n_clust, mc.cores = n_cores, function(i) {
      DefaultAssay(seurat_obj) <- "RNA"
      ident1 <- cell_names[i]
      ident2 <- cell_names[cell_names != cell_names[i]]

      print(paste(ident1, "vs", paste(ident2, collapse = " ")))
      
      table <- FindMarkers(seurat_obj,
        ident.1 = ident1, ident.2 = ident2, only.pos = TRUE)
      table$Gene.name.uniq <- rownames(table)
      table$cluster <- rep(cell_names[i], nrow(table))
      return(table)
    })

  if(save_raw) {
    save_name <- paste0(file_prefix,
      "_clusters_markers_raw_list_", script_name,"_.RDS")
    saveRDS(marker_results, dataPath(save_name))
    print(paste0("raw table saved at ", dataPath(save_name)))
  }

  marker_results <- dplyr::bind_rows(marker_results)
  marker_subset <- marker_results[marker_results$p_val < pval_cutoff,]

  marker_subset$pete_score <- marker_subset$pct.1 *
    marker_subset$avg_logFC * (marker_subset$pct.1 / marker_subset$pct.2)
  marker_subset <- marker_subset[(order(marker_subset$cluster,
    -1 * (marker_subset$pete_score))),]

  marker_table <- dplyr::inner_join(
    marker_subset, gene_info, by = "Gene.name.uniq")

  if(save_table) {
    save_name <- paste0(file_prefix,
      "_cluster_marker_table_", script_name,"_.RDS")
    saveRDS(marker_table, dataPath(save_name))
    print(paste0("Marker table saved at ", dataPath(save_name)))
  }
  
  print("table dimensions:")
  print(dim(marker_table))
  print("markers per cluster:")
  print(table(marker_table$cluster))
  
  return(marker_table)
}

mcPlotMarkers <- function (seurat_obj, diff_results,
  folder_prefix = "cluster-", n_genes = 100, n_cores = 10) {
  
  if (DefaultAssay(seurat_obj) != "RNA") {
    stop("Default assay is not RNA")
  }
  
  dir.create(figurePath(paste0(folder_prefix, "top-markers/")),
    showWarnings = FALSE)
  
  cell_names <- as.character(unique(Idents(seurat_obj)))
  n_clust <- seq_along(unique(Idents(seurat_obj)))
  
  print("active identities:")
  print(cell_names)

  gene_table <- table(diff_results$cluster)
  seq_nums <- seq(1, n_genes, by = 20)

  parallel::mclapply(n_clust, mc.cores = n_cores,
    function(i) {
      genes <- diff_results[diff_results$cluster == cell_names[i],
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
            folder_prefix, "top-markers/", cell_names[i], "_top_",
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