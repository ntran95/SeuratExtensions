
mcPlotMarkers <- function (marker_table) {
  dir.create(figurePath("top-cluster-markers/"))
  n_genes <- 100

  gene_table <- table(marker_table$cluster)
  seq_nums <- seq(1, n_genes, by = 20)
  cell_names <- levels(marker_table$cluster)

  GenerateFPlotPanels <- function(i) {
    genes <- marker_table[marker_table$cluster == cell_names[i],
      "Gene.name.uniq"]
    for(j in 1:length(seq_nums)) {
        to_plot <- genes[seq_nums[j]:(seq_nums[j] + 19)]
        
        if (NA %in% to_plot) {break}

        f <- FeaturePlot(obj_integrated, to_plot,
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
  }

  n_clusters <- 1:length(gene_table)
  # n_clusters <- c(7:10,12) # exclude non specific cell types
  parallel::mclapply(n_clusters, GenerateFPlotPanels, mc.cores = 10)
  # system("pkill -u ddiaz")
}