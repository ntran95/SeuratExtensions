diffConditionMrkrs <- function(
  seurat_obj, group_clusters = NULL, cell_specific = FALSE,
  n_cores = 1, save_raw = TRUE, pval_cutoff = 0.05, save_collapsed = TRUE,
  file_prefix_list = gsub(":|\ ", "-", Sys.time()),
  file_prefix_collapsed = gsub(":|\ ", "-", Sys.time())) {
  
  if (DefaultAssay(seurat_obj) != "RNA") {
    print("Changing default assay to RNA")
    DefaultAssay(seurat_obj) <- "RNA"
  }

  if (cell_specific) {
    marker_table <- readRDS(dataPath(paste0(
      "Clusters_anchored_DF", script_name,"_.RDS")))
    folder <- paste0("cell_specific-diff-", gsub(":|\ ", "-", Sys.time()),"/")
  } else {
    folder <- paste0("cell_type-diff-", gsub(":|\ ", "-", Sys.time()),"/")
  }

  dir.create(figurePath(folder), showWarnings = FALSE)

  all_BCs <- Idents(seurat_obj)
  cell_type <- as.character(unique(Idents(seurat_obj)))
  meta <- seurat_obj@meta.data

  trt <- unique(meta$data.set)
  trt_cnt <- seq(1, (length(trt) * 2) - 1, by = 2)
  trt_BCs <- rep(list(character()), length(trt) * 2)

  if (is.null(group_clusters)) {
    n_clust <- seq_along(unique(Idents(seurat_obj)))
    iterations <- n_clust
  } else {
    n_clust <- which(unique(Idents(seurat_obj)) %in% group_clusters)
    iterations <- 1
  }

  marker_list <- list()[n_clust]
  table_list <- parallel::mclapply(iterations, mc.cores = n_cores,
    FUN = function(h) {
      print(paste0("Number of apply iterations = ",
        iterations, " over ", n_cores, " cores"))
      print(cell_type[n_clust])
      askYesNo("Proceed?")

      if (is.null(group_clusters)) {
        cell_group <- cell_type[h]
      } else {
        cell_group <- cell_type[n_clust]
      }

      for(i in seq_along(trt_cnt)) {
        trt_BCs[[trt_cnt[i]]] <- rownames(meta[
          meta$data.set == trt[i] & meta$cell.type.ident %in% cell_group,])
        trt_BCs[[trt_cnt[i] + 1]] <- rownames(meta[
          !meta$data.set == trt[i] & meta$cell.type.ident %in% cell_group,])
      }

      diff_results <- list() # initalize list

      if (cell_specific) {
        cell_type_markers <- marker_table$Gene.name.uniq[
          marker_table$cell.type.ident %in% cell_group]
      } else {
        cell_type_markers <- NULL
      }

      for (i in seq_along(trt_cnt)) {
        print("Running FindMarkers:")
        print(paste(trt[i], "vs. other timepoints"))

        condition1 <- trt_BCs[[trt_cnt[i]]]
        condition2 <- trt_BCs[[trt_cnt[i] + 1]]
        
        diff_results[[i]] <- FindMarkers(seurat_obj, only.pos = FALSE,
          ident.1 = condition1, ident.2 = condition2, logfc.threshold = 0.25,
          features = cell_type_markers)
      }

      for (i in seq_along(trt_cnt)) {
        diff_results[[i]] <-
          diff_results[[i]][diff_results[[i]]$p_val < pval_cutoff,]
        
        diff_results[[i]]$pete_score <-
          diff_results[[i]]$pct.1 * diff_results[[i]]$avg_logFC * (
            diff_results[[i]]$pct.1 / diff_results[[i]]$pct.2)

        diff_results[[i]] <- diff_results[[i]][
          order(-1 * (diff_results[[i]]$pete_score)),]

        diff_results[[i]]$cell.type.and.trt <-
          paste0(paste0(cell_group, "_", trt[i]), collapse = "-")
        
        diff_results[[i]]$Gene.name.uniq <- rownames(diff_results[[i]])
        diff_results[[i]] <- inner_join(diff_results[[i]],
          gene_info, by = "Gene.name.uniq")
      }
      return(diff_results)
    }
  ) # end mclapply

  if (cell_specific) {
    list_name <- paste0(file_prefix_list, "_marker_list_specific_")
    collapsed_name <- paste0(file_prefix_collapsed, "_all_markers_specific_")
  } else {
    list_name <- paste0(file_prefix_list, "_marker_list_")
    collapsed_name <- paste0(file_prefix_collapsed, "_all_markers_")
  }

  if (save_raw) {
    saveRDS(table_list, dataPath(paste0(list_name, script_name,"_.RDS")))
    table_list <- readRDS(dataPath(paste0(list_name, script_name,"_.RDS")))
  }

  # Collapse list of lists
  if (is.null(group_clusters)) {
    n_results <- n_clust
  } else {
    n_results <- 1
  }

  all_markers <- lapply(n_results, function (i) {bind_rows(table_list[[i]])})
  all_markers <- bind_rows(all_markers)

  if (save_collapsed) {
    saveRDS(all_markers,
      dataPath(paste0(collapsed_name, script_name,"_.RDS")))
    all_markers <- readRDS(
      dataPath(paste0(collapsed_name, script_name,"_.RDS")))
  }
  
  print("marker list (multiple data frames) path:")
  print(dataPath(paste0(list_name, script_name,"_.RDS")))
  cat("\n")
  print("collapsed marker list (single data frame) path:")
  print(dataPath(paste0(collapsed_name, script_name,"_.RDS")))
  
  return(all_markers)
}


# ==== Plot results from diffConditionClust
diffConditionPlots <- function(seurat_obj, input_file = NULL,
  folder_prefix = gsub(":|\ ", "-", Sys.time()), short_sig_figs = TRUE,
  n_genes = 200, n_cores = 4, split = FALSE) {

  if (is.null(input_file)) {
    stop(paste0("Input file with columns Gene.name.uniq ",
      "'cell.type.and.trt' and 'cell.type.ident' required"))
  }
  
  all_markers <- input_file
  print("Number or results for each combination of treatment and cell type:")
  print(table(all_markers$cell.type.and.trt))

  if (short_sig_figs) {
    seurat_obj@meta.data$cell.type.and.trt <- paste0(
      seurat_obj@meta.data$cell.type.ident, "_",
      seurat_obj@meta.data$data.set)

    # Shortening sig figs for plotting
    num_cols <- which(lapply(all_markers, function(i){class(i)}) == "numeric")
    all_markers[,num_cols] <- signif(all_markers[,num_cols], 3)

    DefaultAssay(seurat_obj) <- "RNA"
    assay <- "assay_RNA_"
  }

  seq_nums <- seq(1, n_genes, by = 20)
  meta <- seurat_obj@meta.data

  # get index change for cell type and treatment
  ind_chng <- match(unique(all_markers$cell.type.and.trt),
    all_markers$cell.type.and.trt)

  if (split){
    dir.create(figurePath(paste0(folder_prefix, "-vln-plots")),
      showWarnings = FALSE)
    dir.create(figurePath(paste0(folder_prefix, "-feat-plots")),
      showWarnings = FALSE)
  } else {
    dir.create(figurePath(paste0(folder_prefix, "-feat-vln-plots")),
      showWarnings = FALSE)
  }

parallel::mclapply(seq_along(ind_chng), mc.cores = n_cores, 
    function (i) {

      ifelse(seq_along(ind_chng[i]) == tail(seq_along(ind_chng),1),
        ind_range <- ind_chng[i]:((ind_chng[i + 1]) - 1),
        ind_range <- ind_chng[i]:nrow(all_markers))
      
      population <- paste0(all_markers$cell.type.and.trt[ind_range], " vs. ",
        strsplit(all_markers$cell.type.and.trt[ind_range],"-")[[1]][1],"-other")
      
      stats <- paste0("p-value: ", all_markers$p_val[ind_range], ", ",
        "avg. logFC: ", all_markers$avg_logFC[ind_range], ", ",
        "pct.1: ", all_markers$pct.1[ind_range], ", ",
        "pct.2: ", all_markers$pct.2[ind_range], ", ",
        "pete_score: ", all_markers$pete_score[ind_range])
      
      genes <- all_markers$Gene.name.uniq[ind_range]
      genes <- genes[1:n_genes]

      for(j in seq_along(seq_nums)) {
        cell_ident <- gsub("_.*","",
          all_markers$cell.type.and.trt[ind_range][1])

        sub_ind <- seq_nums[j]:(seq_nums[j]+19)
        pop_sub <- population[sub_ind]
        stats_sub <- stats[sub_ind]

        trt_with_index <- all_markers$cell.type.and.trt[ind_chng[i]]
        to_plot <- genes[sub_ind]
        if (NA %in% to_plot) {break}

        print(paste0(
          "generating plots ", seq_nums[j],"-", (seq_nums[j]+19)))
        vln_list <- VlnPlot(seurat_obj, to_plot, pt.size = 0.25,
          idents = cell_ident, cols = trt_colors, combine = FALSE,
          group.by = "data.set")
        
        for(k in seq_along(vln_list)) {
          vln_list[[k]] <- vln_list[[k]] + NoLegend() + labs(
            caption = paste(pop_sub[k], "\n", stats_sub[k])) +
          theme(plot.caption = element_text(hjust = 0))
        }

        feat_list <- FeaturePlot(seurat_obj, to_plot, reduction = "umap",
          pt.size = 0.25, combine = FALSE)
        
        for(k in seq_along(feat_list)) {
          feat_list[[k]] <- feat_list[[k]] + NoLegend() + NoAxes() + labs(
            caption = paste(pop_sub[k], "\n", stats_sub[k])) +
          theme(plot.caption = element_text(hjust = 0))
        }

        if(split) {
          vln_path <- figurePath(paste0(folder_prefix, "-vln-plots/",
            trt_with_index, "_top_", seq_nums[j],
            "-", (seq_nums[j] + 19),"_features.png"))

          feat_path <- figurePath(paste0(folder_prefix, "-feat-plots/",
            trt_with_index, "_top_", seq_nums[j],
            "-", (seq_nums[j] + 19),"_features.png"))

          png(vln_path, width = 30, height = 25, units = "in", res = 200)
          print(cowplot::plot_grid(plotlist = vln_list))
          dev.off()
          
          png(feat_path, width = 30, height = 25, units = "in", res = 200)
          print(cowplot::plot_grid(plotlist = feat_list))
          dev.off()
        } else {
          combined_path <- figurePath(paste0(folder_prefix, "-feat-vln-plots/",
            trt_with_index, "_top_", seq_nums[j],
            "-", (seq_nums[j] + 19),"_features.png"))

          combined_list <- list()[1:40]
          seqs <- seq(1, 40, by = 5)

          for (x in seq(2, 8, by = 2)) {
            first_row <- (seqs[x-1]):(seqs[x-1]+4)
            next_row <- (seqs[x]):(seqs[x]+4)
            plot_slice <- ((((x/2)*5)-4):((x/2)*5))

            combined_list[first_row] <- feat_list[plot_slice]
            combined_list[next_row] <- vln_list[plot_slice]
          }

          png(combined_path,
            width = 50, height = 60, units = "in", res = 100)
          print(cowplot::plot_grid(plotlist = combined_list, ncol = 5))
          dev.off()
        }
      } # end print plot loop
    }
  ) # end mclappy vln
}


if (FALSE) {
  log_mat <- seurat_obj@assays$RNA@data

  geneMarkerDist <- function(i, group1) {
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
    gene_marker_dist, group1 = TRUE, mc.cores = n_cores)

  stat_sum2 <- parallel::mclapply(1:nrow(all_markers),
    gene_marker_dist, group1 = FALSE, mc.cores = n_cores)

  save(stat_sum1, stat_sum2,
    file = dataPath(paste0("stat_sums", script_name, "_.RData")))
  load(file = dataPath(paste0("stat_sums", script_name, "_.RData")))
}