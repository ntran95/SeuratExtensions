mcFindMarkers <- function(seurat_obj,
  save_raw = TRUE, p_val_cutoff = 0.05, cores = 10) {
  if (DefaultAssay(seurat_obj) != "RNA") {
  stop("Default assay is not RNA")
  }

  cell_names <- unique(Idents(seurat_obj))
  n_clust <- seq_along(unique(Idents(seurat_obj)))

  marker_results <- list()[n_clust]
  marker_results <- parallel::mclapply(
    n_clust, mc.cores = cores, function(i) {
      DefaultAssay(seurat_obj) <- "RNA"
      ident1 <- cell_names[i]
      ident2 <- cell_names[cell_names != cell_names[i]]
      
      table <- FindMarkers(seurat_obj,
        ident.1 = ident1, ident.2 = ident2, only.pos = TRUE)
      table$Gene.name.uniq <- rownames(table)
      table$cluster <- rep(cell_names[i], nrow(table))
      return(table)
    })

  if(save_raw) {
    save_name <- paste0("clusters_markers_raw_list_", script_name,"_.RDS")
    saveRDS(marker_results, dataPath(save_name))
    print(paste0("raw table saved at ", dataPath(save_name)))
  }

  marker_results <- dplyr::bind_rows(marker_results)
  marker_subset <- marker_results[marker_results$p_val < p_val_cutoff,]
  
  marker_subset$pete_score <- marker_subset$pct.1 *
    marker_subset$avg_logFC * (marker_subset$pct.1 / marker_subset$pct.2)
  marker_subset <- marker_subset[(order(marker_subset$cluster,
    -1 * (marker_subset$pete_score))),]

  marker_table <- dplyr::inner_join(
    marker_subset, gene_info, by = "Gene.name.uniq")
  
  print("table dimensions ")
  print(dim(marker_table))
  print("markers per cluster")
  print(table(marker_table$cluster))
  
  return(marker_table)
}