# ================ Diff expression for every treatment with NEW Cell classes
if (FALSE) {
  if (DefaultAssay(obj_integrated) != "RNA") {
    print("Changing default assay to RNA")
    DefaultAssay(obj_integrated) <- "RNA"
  }

  cell_specific <- TRUE
  if (cell_specific) {
    marker_table <- readRDS(dataPath(paste0(
      "Clusters_anchored_DF", script_name,"_.RDS")))
    folder <- "diff_exp_cell_specific/"
  } else {
    folder <- "diff_exp_by_cell_type/"}

  dir.create(figurePath(folder))

  all_BCs <- Idents(obj_integrated)
  cell_type <- as.character(unique(Idents(obj_integrated)))
  meta <- obj_integrated@meta.data

  trt <- unique(meta$data.set)
  trt_cnt <- seq(1, (length(trt) * 2) - 1, by = 2)
  trt_BCs <- rep(list(character()), length(trt) * 2) # init list

  n_clust <- seq_along(unique(Idents(obj_integrated)))
  marker_list <- list()[n_clust] # init list

  MarkerCombos <- function(h) {

    for(i in seq_along(trt_cnt)) {
      trt_BCs[[trt_cnt[i]]] <- rownames(meta[
        meta$data.set == trt[i] & meta$cell.type.ident == cell_type[h],])
      trt_BCs[[trt_cnt[i] + 1]] <- rownames(meta[
        !meta$data.set == trt[i] & meta$cell.type.ident == cell_type[h],])
    }

    diff_results <- list() # initalize list

    if (cell_specific) {
      cell_type_markers <- marker_table$Gene.name.uniq[
        marker_table$cell.type.ident == cell_type[h]]
    } else {cell_type_markers <- NULL}

    for (i in seq_along(trt_cnt)) {
      condition1 <- trt_BCs[[trt_cnt[i]]]
      condition2 <- trt_BCs[[trt_cnt[i] + 1]]
      diff_results[[i]] <- FindMarkers(obj_integrated, only.pos = FALSE,
        ident.1 = condition1, ident.2 = condition2, logfc.threshold = 0.25,
        features = cell_type_markers)
      print(dim(diff_results[[i]]))
    }

    for (i in seq_along(trt_cnt)) {
      diff_results[[i]] <- diff_results[[i]][diff_results[[i]]$p_val < 0.05,]
      diff_results[[i]]$pct.ratio <- 1:nrow(diff_results[[i]])
      
      diff_results[[i]]$pct.ratio <- 
        diff_results[[i]]$pct.1 / diff_results[[i]]$pct.2
      
      diff_results[[i]] <- diff_results[[i]][(
        order(diff_results[[i]]$pct.ratio, decreasing = TRUE)),]

      diff_results[[i]]$cell.type.and.trt <- paste0(cell_type[h], "_", trt[i])
      diff_results[[i]]$Gene.name.uniq <- rownames(diff_results[[i]])
      diff_results[[i]] <- inner_join(diff_results[[i]],
        gene_info, by = "Gene.name.uniq")
    }

    table_list <- list()[seq_along(trt_cnt)]

    for (i in seq_along(trt_cnt)) {
      finished_table <- diff_results[[i]]
      table_list[[i]] <- finished_table
    }
    return(table_list)
  } 

  marker_list <- parallel::mclapply(n_clust, FUN = MarkerCombos, mc.cores = 10)

  cell_specific <- TRUE
  if (cell_specific) {
    file_name <- "MarkerList_specific_"
  } else {file_name <- "MarkerList_"}

  if (FALSE) {
    saveRDS(marker_list, dataPath(paste0(file_name, script_name,"_.RDS")))
    marker_list <- readRDS(dataPath(paste0(file_name, script_name,"_.RDS")))
  }

  # Collapse list into a single data frame
  n_clust <- seq_along(unique(Idents(obj_integrated)))
  all_markers <- lapply(n_clust, function (i) {bind_rows(marker_list[[i]])})
  all_markers <- bind_rows(all_markers)
  class(all_markers)

  if (FALSE) {
    saveRDS(all_markers, dataPath(paste0(file_name, script_name,"_.RDS")))
    all_markers <- readRDS(dataPath(paste0(file_name, script_name,"_.RDS")))
  }
}

# ================ Plot results from above (each cell type vs all timepoints)
if (FALSE) { # place holder for new function
  cell_specific <- TRUE
  if (cell_specific) {
    file_name <- "MarkerList_specific_"
    all_markers <- readRDS(dataPath(paste0(file_name, script_name,"_.RDS")))
    folder <- "diff_exp_cell_specific/"
  } else {
    file_name <- "MarkerList_"
    all_markers <- readRDS(dataPath(paste0(file_name, script_name,"_.RDS")))
    folder <- "diff_exp_by_cell_type/"
  }

  table(all_markers$cell.type.and.trt) < 20
  length(all_markers$cell.type.and.trt)

  short_sig_figs <- TRUE
  if (short_sig_figs) {
    obj_integrated@meta.data$cell.type.and.trt <- paste0(
      obj_integrated@meta.data$cell.type.ident, "_",
      obj_integrated@meta.data$data.set)

    # Shortening sig figs for plotting
    all_markers$avg_logFC <- signif(all_markers$avg_logFC, 3)
    all_markers$p_val <- signif(all_markers$p_val, 3)
    all_markers$pct.ratio <- signif(all_markers$pct.ratio, 2)
    all_markers$pct.1 <- signif(all_markers$pct.1, 3)
    all_markers$pct.2 <- signif(all_markers$pct.2, 3)

    DefaultAssay(obj_integrated) <- "RNA"
    assay <- "assay_RNA_"
  }
}


# ---- Get summary stats for all diff genes  ---- #
if (FALSE) { # place holder for new function
  summary_stat <- FALSE
  if (summary_stat) {
    log_mat <- obj_integrated@assays$RNA@data

    gene_marker_dist <- function(i, group1) {
      if (group1 == TRUE) {bc_ind <- grepl(
        all_markers$cell.type.and.trt[i], meta$cell.type.and.trt)
        group <- "group-1"
      } else {bc_ind <- !grepl(
        all_markers$cell.type.and.trt[i], meta$cell.type.and.trt)
        group <- "group-2"}
      
      gene_vec <- log_mat[all_markers$Gene.name.uniq[i], bc_ind]
      
      stats <- summary(log_mat[all_markers$Gene.name.uniq[i], bc_ind])
      stats <- c(paste(all_markers$Gene.name.uniq[i], group),
        paste(names(stats), signif(unname(stats), 2), sep = ":"))
      return(stats)
    }
    
    # This takes awhile to run for the macrophage data, even with 10 cores
    stat_sum1 <- parallel::mclapply(1:nrow(all_markers),
      gene_marker_dist, group1 = TRUE, mc.cores = 10)

    stat_sum2 <- parallel::mclapply(1:nrow(all_markers),
      gene_marker_dist, group1 = FALSE, mc.cores = 10)

    save(stat_sum1, stat_sum2,
      file = dataPath(paste0("stat_sums", script_name, "_.RData")))
    load(file = dataPath(paste0("stat_sums", script_name, "_.RData")))
  }

  # get index change for cell type and treatment
  ind_chng <- match(unique(all_markers$cell.type.and.trt),
    all_markers$cell.type.and.trt)


  # ---- Violin Plot panels ---- #
  dir.create(figurePath(paste0(folder, "Vplot-top-markers-by-treatment/")))
  n_genes <- 200
  seq_nums <- seq(1, n_genes, by = 20)
  meta <- obj_integrated@meta.data

  parallel::mclapply(seq_along(ind_chng), mc.cores = 10, function (i) {
    ifelse(seq_along(ind_chng[i]) == tail(seq_along(ind_chng),1),
      ind_range <- ind_chng[i]:((ind_chng[i + 1]) - 1),
      ind_range <- ind_chng[i]:nrow(all_markers))
    
    population <- paste0(all_markers$cell.type.and.trt[ind_range], " vs. ",
      (strsplit(all_markers$cell.type.and.trt[ind_range],"-")[[1]][1]), "-other")
    
    stats <- paste0("p-value: ", all_markers$p_val[ind_range], ", ",
      "avg. logFC: ", all_markers$avg_logFC[ind_range], ", ",
      "pct.1: ", all_markers$pct.1[ind_range], ", ",
      "pct.2: ", all_markers$pct.2[ind_range], ", ",
      "pct.ratio: ", all_markers$pct.ratio[ind_range])
    
    genes <- all_markers$Gene.name.uniq[ind_range]
    genes <- genes[1:n_genes]

    for(j in seq_along(seq_nums)) {
      cell_ident <- gsub("_.*","", all_markers$cell.type.and.trt[ind_range][1])

      sub_ind <- seq_nums[j]:(seq_nums[j]+19)
      pop_sub <- population[sub_ind]
      stats_sub <- stats[sub_ind]

      to_plot <- genes[sub_ind]
      if (NA %in% to_plot) {break}

      v <- VlnPlot(obj_integrated, to_plot, pt.size = 0.25,
        idents = cell_ident, cols = trt_colors,
        combine = FALSE, group.by = "data.set")
      for(k in seq_along(v)) {
        v[[k]] <- v[[k]] + NoLegend() + labs(
          caption = paste(pop_sub[k], "\n", stats_sub[k])) +
        theme(plot.caption = element_text(hjust = 0))
      }

      path <- figurePath(paste0(folder, "Vplot-top-markers-by-treatment/",
        all_markers$cell.type.and.trt[ind_chng[i]], "_top_", seq_nums[j],
        "-", (seq_nums[j] + 19),"_features.png"))

      png(path, width = 30, height = 25, units = "in", res = 200)
      print(cowplot::plot_grid(plotlist = v))
      dev.off()
    }
  })


  # ---- Feature Plot panels ---- #
  dir.create(figurePath(paste0(folder, "Fplot-top-markers-treatment/")))
  n_genes <- 200
  seq_nums <- seq(1, n_genes, by = 20)

  parallel::mclapply(seq_along(ind_chng), mc.cores = 12, function (i) {
    ifelse(seq_along(ind_chng[i]) == tail(seq_along(ind_chng),1),
      ind_range <- ind_chng[i]:((ind_chng[i + 1]) - 1),
      ind_range <- ind_chng[i]:nrow(all_markers))
    
    population <- paste0(all_markers$cell.type.and.trt[ind_range], " vs. ",
      (strsplit(all_markers$cell.type.and.trt[ind_range],"-")[[1]][1]), "-other")
    
    stats <- paste0("p-val: ", all_markers$p_val[ind_range], ", ",
      "avg. logFC: ", all_markers$avg_logFC[ind_range], ", ",
      "pct.1: ", all_markers$pct.1[ind_range], ", ",
      "pct.2: ", all_markers$pct.2[ind_range], ", ",
      "pct.ratio: ", all_markers$pct.ratio[ind_range])

    genes <- all_markers$Gene.name.uniq[ind_range]
    genes <- genes[1:n_genes]

    for(j in seq_along(seq_nums)) {
      sub_ind <- seq_nums[j]:(seq_nums[j] + 19)
      to_plot <- genes[sub_ind]
      pop_sub <- population[sub_ind]
      stats_sub <- stats[sub_ind]
      
      if (NA %in% to_plot) {break}
      f <- FeaturePlot(obj_integrated, to_plot,
        reduction = "umap", pt.size = 0.25, combine = FALSE)
      
      for(k in seq_along(f)) {
        f[[k]] <- f[[k]] + NoLegend() + NoAxes() + labs(
          caption = paste(pop_sub[k], "\n", stats_sub[k])) +
        theme(plot.caption = element_text(hjust = 0))
      }

      path <- figurePath(paste0(folder, "Fplot-top-markers-treatment/",
        all_markers$cell.type.and.trt[ind_chng[i]], "_top_", seq_nums[j],
        "-", (seq_nums[j] + 19),"_features.png"))
      png(path, width = 30, height = 25, units = "in", res = 200)
      print(cowplot::plot_grid(plotlist = f))
      dev.off()
    }
  })
}
