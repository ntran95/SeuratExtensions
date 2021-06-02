if (FALSE) { 
# ================ Diff Expression on anchored data
DefaultAssay(obj_integrated) <- "RNA"

# Parallel solution for FindAllMarkers
cell_names <- unique(Idents(obj_integrated))
n_clust <- seq_along(unique(Idents(obj_integrated)))

mcFindMarkers <- function(i) {
  ident1 <- cell_names[i]
  ident2 <- cell_names[cell_names != cell_names[i]]
  table <- FindMarkers(obj_integrated,
    ident.1 = ident1, ident.2 = ident2, only.pos = TRUE)
  table$Gene.name.uniq <- rownames(table)
  table$cell.type.ident <- rep(cell_names[i], nrow(table))
  return(table)
}

marker_results <- list()[n_clust]
marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = 10)

if(TRUE) {
  saveRDS(marker_results, dataPath(paste0(
    "Clusters_anchored_cell_type_", script_name,"_.RDS")))
  markers <- readRDS(dataPath(paste0(
    "Clusters_anchored_cell_type_", script_name,"_.RDS")))
}

markers <- dplyr::bind_rows(markers) # turn into a single data frame

marker_subset <- markers[markers$p_val < 0.05,]
marker_subset$pct.ratio <- 1:nrow(marker_subset)
marker_subset$pct.ratio <- marker_subset$pct.1 / marker_subset$pct.2
marker_subset <- marker_subset[(order(marker_subset$cell.type.ident,
  -1 * (marker_subset$pct.ratio))),]

dim(marker_subset)
table(marker_subset$cell.type.ident)
marker_table <- inner_join(marker_subset, gene_info, by = "Gene.name.uniq")

WriteXLS::WriteXLS(marker_table, figurePath(paste0("cluster_markers_tree_ident",
  assay,".xlsx")), row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)


# ================ mclapply for plotting
all_BCs <- Idents(obj_integrated)
cell_type <- as.character(unique(Idents(obj_integrated)))
n_clust <- seq_along(unique(Idents(obj_integrated)))

dir.create(figurePath(folder))
n_genes <- 20
gene_table <- table(marker_table$cluster)
seq_nums <- seq(1, n_genes, by = 20)

GenerateFPlotPanels <- function(i) {
  genes <- marker_table[
    marker_table$cluster == cell_type[i], "Gene.name.uniq"]
  for(j in 1:length(seq_nums)) {
      to_plot <- genes[seq_nums[j]:(seq_nums[j] + 19)]
      
      if (NA %in% to_plot) {break}

      f <- FeaturePlot(obj_integrated, to_plot, slot = "data",
        reduction = "umap", pt.size = 0.25, combine = FALSE)
      for(k in 1:length(f)) {
        f[[k]] <- f[[k]] + NoLegend() + NoAxes()
      }

      path <- figurePath(paste0(folder, "_", cell_type[i],
        "_top_", seq_nums[j], "-", (seq_nums[j] + 19),"_features.png"))
      png(path, width = 30, height = 25, units = "in", res = 200)
      print(cowplot::plot_grid(plotlist = f))
      dev.off()
  }
}

n_clusters <- 1:length(gene_table)
ptm <- proc.time()
parallel::mclapply(n_clusters, GenerateFPlotPanels, mc.cores = 10)
time_diff <- proc.time() - ptm
time_diff
}