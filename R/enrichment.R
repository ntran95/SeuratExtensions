# read in single marker table
enrichment_analysis <- function(marker_tbl,
                                dir_name,
                                minGSSize = 3,
                                maxGSSize = 800,
                                organism = "org.Dr.eg.db",
                                reactome_organism  = "zebrafish",
                                kegg_organism  = "dre",
                                pvalueCutoff = 0.05,
                                showCategory = 15,
                                lower_bound_logFC = -1.0,
                                upper_bound_logFC = 1.0){
  #if ident column is empty, call on modify_marker_tbl(reprocess = TRUE)
  if ("ident" %notin% colnames(marker_tbl)){
    print("reprocessing marker table")
    
    marker_tbl <- modify_marker_tbl(marker_tbl,
                                    findAllMarkers = FALSE,
                                    ident.1 = unique(marker_tbl$cluster),
                                    ident.2 = NULL,
                                    pval_cutoff = 0.05,
                                    reprocess = TRUE)
  }
  
  # create a sub directory based on ident column of marker_tbl
  dir_name_sub <- unique(marker_tbl$ident)
  
  dir.create(figurePath(paste0(dir_name, "/", dir_name_sub)),
             showWarnings = F,recursive = T)
  
  # make go obj
  gse <- SeuratExtensions::make_gsea_obj(marker_tbl = marker_tbl,
                                         minGSSize = minGSSize,
                                         maxGSSize = maxGSSize,
                                         organism = organism,
                                         pvalueCutoff = pvalueCutoff)
  
  # make kegg obj
  kegg <- SeuratExtensions::make_kegg_obj(marker_tbl = marker_tbl,
                                          minGSSize = minGSSize,
                                          maxGSSize = maxGSSize,
                                          organism = organism,
                                          kegg_organism  = "dre",
                                          pvalueCutoff = pvalueCutoff)
  
  # make reactome obj
  #reactome pathway enrichment
  reactome <-
    SeuratExtensions::make_reactome_obj(marker_tbl = marker_tbl,
                                        organism = organism,
                                        reactome_organism  = reactome_organism,
                                        pvalueCutoff = pvalueCutoff,
                                        minGSSize = minGSSize)
  
  # append go terms to marker tbl
  marker_tbl_sub_gse <-
    SeuratExtensions::add_gsea_to_marker_tbl(marker_tbl = marker_tbl,
                                             gsea = gse)
  
  # append kegg terms to marker tbl
  
  marker_tbl_sub_kegg <-
    add_kegg_to_marker_tbl(marker_tbl = marker_tbl,
                           kegg = kegg)
  
  # append reactome terms to marker tbl
  marker_tbl_reactome <-
    SeuratExtensions::add_reactome_to_marker_tbl(marker_tbl = marker_tbl,
                                                 reactome = reactome)
  
  # merge go and kegg columms to single marker tbl
  marker_tbl_gse_kegg <- inner_join(marker_tbl_sub_gse,
                                    marker_tbl_sub_kegg)
  
  # merge go kegg reactome columsn to single marker tbl, complete_gsea_marker_tbl
  complete_gsea_marker_tbl <- inner_join(marker_tbl_gse_kegg,
                                         marker_tbl_reactome)
  
  # store obj in list
  res_list <- list(gse,kegg,reactome)
  names(res_list) <- c("GSEA", "KEGG", "REACTOME")
  
  # iterate through obj list, graph enrich plots
  #check that results show > 0 enriched terms before printing to file
  for (j in seq_along(res_list)){
    if (nrow(res_list[[j]]@result) != 0 ){
      #make plots, Reactome
      enrich_plots <-
        SeuratExtensions::make_enrich_plot(enrich_obj = res_list[[j]],
                                           marker_tbl = complete_gsea_marker_tbl,
                                           is.enrich.reactome = F,
                                           analysis.type = names(res_list)[j],
                                           showCategory = showCategory,
                                           dir_name = dir_name,
                                           dir_name_sub = dir_name_sub,
                                           save = T,
                                           x_axis = "GeneRatio")
      
      
      
    }else{
      next
    } 
  }
  
  
  
  # graph volcano plots
  #create volcano plot
  vp <-
    SeuratExtensions::volcano_plot(marker_tbl = marker_tbl,
                                   ident = unique(marker_tbl$ident),
                                   gene_col_name = "Gene.name.uniq",
                                   lower_bound_logFC = -1.0,
                                   upper_bound_logFC = 1.0,
                                   logFC_hline = c(-1.0, 1.0)
    )
  
  
  # make interactive vp
  vp_label_all <-
    SeuratExtensions::volcano_plot(marker_tbl = marker_tbl,
                                   ident = unique(marker_tbl$ident),
                                   gene_col_name = "Gene.name.uniq",
                                   lower_bound_logFC = -1.0,
                                   upper_bound_logFC = 1.0,
                                   logFC_hline = c(-1.0, 1.0),
                                   label.all = TRUE)
  
  vp.plotly <- plotly::ggplotly(vp_label_all)
  
  htmlwidgets::saveWidget(vp.plotly,
                          figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                                            "interactive_volcano_plot_",
                                            dir_name_sub,
                                            ".html")))
  
  # # featureplots
  # n_genes <- 50
  # seq_nums <- seq(1, n_genes, by = 20)
  # n_ident <- seq_along(unique(marker_tbl$ident))
  # folder_prefix <- "feature-plot"
  # 
  # dir.create(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
  #                              folder_prefix, "-top-markers/")),
  #            showWarnings = FALSE)
  # for (i in seq_along(n_ident)) {
  #   #order by pete_score in ascending order 
  #   #marker_tbl <-marker_tbl[order(marker_tbl$pete_score,decreasing = T),]
  #   
  #   genes <- marker_tbl$Gene.name.uniq
  #   
  #   stats <- paste0("p-value: ", signif(marker_tbl$p_val,digits = 3), ", ",
  #                   "avg. logFC: ", signif(marker_tbl$avg_logFC,digits = 3), ", ",
  #                   "pct.1: ", marker_tbl$pct.1, ", ",
  #                   "pct.2: ", marker_tbl$pct.2, ", ",
  #                   "pete_score: ", signif(marker_tbl$pete_score,3))
  #   population <- marker_tbl$ident
  #   
  #   for(j in 1:length(seq_nums)) {
  #     sub_ind <- seq_nums[j]:(seq_nums[j]+19)
  #     pop_sub <- population[sub_ind]
  #     stats_sub <- stats[sub_ind]
  #     
  #     to_plot <- genes[seq_nums[j]:(seq_nums[j] + 19)]
  #     
  #     if (NA %in% to_plot) {break}
  #     
  #     f <- Seurat::FeaturePlot(obj_integrated, to_plot,
  #                              reduction = "umap", pt.size = 0.25, combine = FALSE)
  #     for(k in 1:length(f)) {
  #       f[[k]] <- f[[k]] + NoLegend() + NoAxes() +
  #         labs( caption = paste(pop_sub[k], "\n", stats_sub[k])) +
  #         theme(plot.caption = element_text(hjust = 0))
  #     }
  #     
  #     path <- figurePath(paste0(dir_name, "/", dir_name_sub, "/",
  #                               folder_prefix, "-top-markers/", 
  #                               unique(marker_tbl$ident), "_top_",
  #                               seq_nums[j], "-", (seq_nums[j] + 19), 
  #                               "_features.png"))
  #     
  #     png(path, width = 40, height = 35, units = "in", res = 300)
  #     print(cowplot::plot_grid(plotlist = f))
  #     dev.off()
  #   } 
  # }
  # 
  # write out complete_gsea_marker_tbl to excel
  openxlsx::write.xlsx(complete_gsea_marker_tbl,
                       file = figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                                                "gene_list_", dir_name_sub, ".xlsx")))


  # save obj list
  save(list = c("gse", "kegg", "reactome"),
       file = dataPath(paste0("saved_gse_obj_",
                              dir_name_sub, "_",
                              script_name, ".RData")))


  
}

