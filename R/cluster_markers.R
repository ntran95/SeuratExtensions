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
  folder_prefix = "cluster-", n_genes = 100, n_cores = 10, 
  split.by = NULL,
  split.by_folder_name = NULL) {
  
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
          
          if(!is.null(split.by)){
            dir.create(figurePath(paste0(folder_prefix, "top-markers/", 
                                         split.by_folder_name, "/")))
            
            dir.create(figurePath(paste0(folder_prefix, "top-markers/", 
                                         split.by_folder_name, "/violin_plots/")),
                       showWarnings = FALSE)
            
            dir.create(figurePath(paste0(folder_prefix, "top-markers/", 
                                         split.by_folder_name, "/feature_plots/")),
                       showWarnings = FALSE)
            
            v  <- VlnPlot(obj_integrated, features = to_plot, 
                            split.by = "data.set", split.plot = TRUE,
                            combine = FALSE)
            
            for(k in 1:length(f)) {
              v[[k]] <- v[[k]] 
            }
            
            v_path <- figurePath(paste0(
              folder_prefix, "top-markers/", split.by_folder_name, 
              "/violin_plots/", 
              cell_names[i], "_top_",
              seq_nums[j], "-", (seq_nums[j] + 19), "_violin.png"))
            
            png(v_path, width = 30, height = 25, units = "in", res = 200)
            print(cowplot::plot_grid(plotlist = v))
            dev.off()  
            
            f <- FeaturePlot(obj_integrated, to_plot,
                             reduction = "umap", pt.size = 0.25, 
                             split.by = "data.set") 
           
            f_path <- figurePath(paste0(
              folder_prefix, "top-markers/", split.by_folder_name, 
              "/feature_plots/", 
              cell_names[i], "_top_",
              seq_nums[j], "-", (seq_nums[j] + 19), "_features.png"))
            
            png(f_path, width = 30, height = 70, units = "in", res = 200)
            print(f)
            dev.off() 
            
            }
      }
      return(cell_names[i])
      print("Done")
    }
  )
}

# near same plotting function as mcPlotMarkers() that does not use mclappy()
top_n_featureplots <- function(seurat_obj,
                              marker_tbl, 
                              n_genes = 50, 
                              dir_name = NULL,
                              folder_prefix = "feature-plot",
                              shape.by = NULL,
                              label = TRUE,
                              group.by = "tree.ident",
                              split.by = NULL){
  # featureplots
  n_genes <- n_genes
  seq_nums <- seq(1, n_genes, by = 20)
  n_ident <- seq_along(unique(marker_tbl$ident))

  #group cells by idents
  Idents(seurat_obj) <- group.by

  if(!is.null(dir_name)){
    dir.create(figurePath(paste0(dir_name)),
             showWarnings = FALSE)
  }

  dir.create(figurePath(paste0(dir_name, "/",
                               folder_prefix, "-top-markers/")),
             showWarnings = FALSE)
  for (i in seq_along(n_ident)) {
    #order by pete_score in ascending order 
    #marker_tbl <-marker_tbl[order(marker_tbl$pete_score,decreasing = T),]
    
    genes <- marker_tbl$Gene.name.uniq
    
    stats <- paste0("p-value: ", signif(marker_tbl$p_val,digits = 3), ", ",
                    "avg. logFC: ", signif(marker_tbl$avg_logFC,digits = 3), ", ",
                    "pct.1: ", marker_tbl$pct.1, ", ",
                    "pct.2: ", marker_tbl$pct.2, ", ",
                    "pete_score: ", signif(marker_tbl$pete_score,3))
    population <- marker_tbl$ident
    
    for(j in 1:length(seq_nums)) {
      sub_ind <- seq_nums[j]:(seq_nums[j]+19)
      pop_sub <- population[sub_ind]
      stats_sub <- stats[sub_ind]
      
      to_plot <- genes[seq_nums[j]:(seq_nums[j] + 19)]
      
      if (NA %in% to_plot) {break}
      
      f <- Seurat::FeaturePlot(seurat_obj, to_plot,
                               reduction = "umap", 
                               pt.size = 0.25, 
                               combine = FALSE,
                               split.by = split.by,
                               label = label)
      for(k in 1:length(f)) {
        f[[k]] <- f[[k]] + NoLegend() + NoAxes() +
          labs( caption = paste(pop_sub[k], "\n", stats_sub[k])) +
          theme(plot.caption = element_text(hjust = 0))
      }
      
      path <- figurePath(paste0(dir_name, "/", folder_prefix, "-top-markers/", 
                                unique(marker_tbl$ident), "_top_",
                                seq_nums[j], "-", (seq_nums[j] + 19), 
                                "_features.png"))
      
      png(path, width = 40, height = 35, units = "in", res = 300)
      print(cowplot::plot_grid(plotlist = f))
      dev.off()
    } 
  }
}

top_n_dotplot <- function(seurat_obj,
                              marker_tbl, 
                              n_genes = 50, 
                              dir_name = NULL,
                              folder_prefix = "dot-plot",
                              shape.by = NULL,
                              label = TRUE,
                              group.by = "tree.ident",
                              split.by = NULL){
  n_genes <- n_genes
  seq_nums <- seq(1, n_genes, by = 20)
  n_ident <- seq_along(unique(marker_tbl$ident))

  #group cells by idents
  Idents(seurat_obj) <- group.by

  if(!is.null(dir_name)){
    dir.create(figurePath(paste0(dir_name)),
             showWarnings = FALSE)
  }

  dir.create(figurePath(paste0(dir_name, "/",
                               folder_prefix, "-top-markers/")),
             showWarnings = FALSE)
  for (i in seq_along(n_ident)) {
    #order by pete_score in ascending order 
    #marker_tbl <-marker_tbl[order(marker_tbl$pete_score,decreasing = T),]
    
    genes <- marker_tbl$Gene.name.uniq
    
    stats <- paste0("p-value: ", signif(marker_tbl$p_val,digits = 3), ", ",
                    "avg. logFC: ", signif(marker_tbl$avg_logFC,digits = 3), ", ",
                    "pct.1: ", marker_tbl$pct.1, ", ",
                    "pct.2: ", marker_tbl$pct.2, ", ",
                    "pete_score: ", signif(marker_tbl$pete_score,3))
    population <- marker_tbl$ident
    
    for(j in 1:length(seq_nums)) {
      sub_ind <- seq_nums[j]:(seq_nums[j]+19)
      pop_sub <- population[sub_ind]
      stats_sub <- stats[sub_ind]
      
      to_plot <- genes[seq_nums[j]:(seq_nums[j] + 19)]
      
      if (NA %in% to_plot) {break}
      
      f <- Seurat::DotPlot(seurat_obj,
                           features = rev(to_plot), #or else genes start at bottom
                           group.by = group.by,
                           split.by = split.by,
                           cols = "RdYlBu")

      f <- f + coord_flip() + theme(
        axis.text.x = element_text(angle = 90, hjust = 1))
      # for(k in 1:length(f)) {
      #   f[[k]] <- f[[k]] + 
      #     labs( caption = paste(pop_sub[k], "\n", stats_sub[k])) +
      #     theme(plot.caption = element_text(hjust = 0))
      # }
      
      path <- figurePath(paste0(dir_name, "/", folder_prefix, "-top-markers/", 
                                unique(marker_tbl$ident), "_top_",
                                seq_nums[j], "-", (seq_nums[j] + 19), 
                                "_dotplot.png"))
      
      png(path, width = 12, height = 10, units = "in", res = 300)
      print(f)
      dev.off()
    } 
  }
}

modify_marker_tbl <- function(marker_tbl, 
                              findAllMarkers = FALSE, 
                              ident.1,
                              ident.2 = NULL,
                              pval_cutoff = 0.05, 
                              reprocess = FALSE){
  if(reprocess){
    #calculate pete_score
    marker_tbl$pete_score <-
      marker_tbl$pct.1 * marker_tbl$avg_logFC * (
        marker_tbl$pct.1 / marker_tbl$pct.2)
    
    #pvalue cutoff
    marker_tbl <- marker_tbl[marker_tbl$p_val < pval_cutoff,]
    #specify which groups to compare to
    if (is.null(ident.2) == TRUE){
      marker_tbl$ident <- paste0(ident.1, "_vs_", "all_other_cells")
    }else{
      marker_tbl$ident <- paste0(ident.1, "_vs_", ident.2)
    }
  
  }else{
    if(findAllMarkers ==TRUE){
      marker_tbl <- merge(marker_tbl, gene_info,
                          by.x = "gene",
                          by.y ="Gene.name.uniq", 
                          all.x = TRUE)
      marker_tbl <- marker_tbl[order(marker_tbl$cluster, marker_tbl$avg_logFC,
                                     decreasing = TRUE),]
      
      marker_tbl <- marker_tbl[order(marker_tbl$cluster), ]
      
      
      names(marker_tbl)[names(marker_tbl) == "gene"] <- "Gene.name.uniq"
      
    }else{
      marker_tbl$Gene.name.uniq <- rownames(marker_tbl)
      marker_tbl <- merge(marker_tbl, gene_info,
                          by.x = "Gene.name.uniq",
                          by.y ="Gene.name.uniq", 
                          all.x = TRUE)
      marker_tbl <- marker_tbl[order(marker_tbl$avg_logFC,
                                     decreasing = TRUE),]
    }
    
    #calculate pete_score
    marker_tbl$pete_score <-
      marker_tbl$pct.1 * marker_tbl$avg_logFC * (
        marker_tbl$pct.1 / marker_tbl$pct.2)
    
    #pvalue cutoff
    marker_tbl <- marker_tbl[marker_tbl$p_val < pval_cutoff,]
    
    
    # #if gene symbol column == ensembl ID, append ensembl id to Gene.stable.ID as well
    # marker_tbl$Gene.stable.ID <- ifelse(is.na(marker_tbl$Gene.stable.ID) == TRUE, 
    #                                     as.character(marker_tbl$Gene.name.uniq),  #if == True
    #                                     as.character(marker_tbl$Gene.stable.ID)) #keep same 
    # 
    # #if gene symbol column == ensembl ID, append ensembl id to Gene.stable.ID as well
    # marker_tbl$Gene.name <- ifelse(is.na(marker_tbl$Gene.name) == TRUE, 
    #                                     as.character(marker_tbl$Gene.name.uniq),  #if == True
    #                                     as.character(marker_tbl$Gene.stable.ID)) #keep same 
    
    
    #specify which groups to compare to
    if (findAllMarkers){
      marker_tbl$ident <- 
        paste0(marker_tbl$cluster, "_vs_", "all_other_cells")
    }else if(is.null(ident.2)== TRUE){
      marker_tbl$ident <- paste0(ident.1, "_vs_", "all_other_cells")
    }else{
      marker_tbl$ident <- paste0(ident.1, "_vs_", ident.2)
    }
  }
  return(marker_tbl)
}