hmap <- function(seurat_obj, genes, 
                 group.by = "cell.type.ident"){
  dotplot <- Seurat::DotPlot(seurat_obj, features = genes,
                             group.by = group.by)
  
  plot_df <- dotplot$data
  
  
  g <-  ggplot(plot_df) +
    geom_tile(aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1),
              color = "gray", size = 1) +
    scale_fill_distiller(
      palette = "RdYlBu") +
    hrbrthemes::theme_ipsum(base_family = "sans", grid = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
          axis.title.y.right = element_text(size=13),
          panel.spacing = unit(.35, "lines"),
          strip.text.x  = element_text(vjust = 0.5, hjust=.5,size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    ylim(rev(levels(plot_df$features.plot))) 
  
  if(group.by == "cell.type.ident.by.data.set"){
    plot_df$groupIdent <- gsub("(.+?)(\\_.*)", "\\1",plot_df$id)
    plot_df$groupIdent <- factor(plot_df$groupIdent,
                                 levels=levels(seurat_obj$cell.type.ident))
    
    g <-  g  %+% plot_df +
      facet_grid( ~ groupIdent, scales='free_x') 
    
  }else if(group.by == "cell.type.and.trt"){
    plot_df$groupIdent <- gsub("(.+?)(\\..*)", "\\2",plot_df$id)
    plot_df$groupIdent <- gsub("^.", "",  plot_df$groupIdent)
    plot_df$groupIdent <- factor(plot_df$groupIdent,
                                 levels=levels(seurat_obj$cell.type.ident))
    g <-  g  %+% plot_df +
      facet_grid( ~ groupIdent, scales='free_x') 
  } else if(group.by == "effector.and.trt"){
    plot_df$groupIdent <- gsub("(.+?)(\\..*)", "\\2",plot_df$id)
    plot_df$groupIdent <- gsub("^.", "",  plot_df$groupIdent)
    plot_df$groupIdent <- factor(plot_df$groupIdent,
                                levels=levels(seurat_obj$effector_groups))
    g <-  g  %+% plot_df +
      facet_grid( ~ groupIdent, scales='free_x') 
  }
  # if(facet_title == "data.set"){
  #   
  #   plot_df$groupIdent <- gsub(".*_", "",plot_df$id)
  #   plot_df$groupIdent <- factor(plot_df$groupIdent,
  #                                levels=levels(seurat_obj$data.set))
  #   
  #   g <-  g  %+% plot_df +
  #     facet_grid( ~ groupIdent, scales='free_x') 
  # }
  
  
  return(g)
  
}

indv_cell_hmap <- function(seurat_obj, genes, group.by = "cell.type.ident"){
  col.min = -2.5
  col.max = 2.5
  
  cells <- colnames(x = seurat_obj)
  
  data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
    object = seurat_obj, slot = "data")[genes, cells, drop = FALSE])))
  
  
  data <- scale(data)
  data <- as.data.frame(MinMax(data = data, min = col.min, max = col.max))
  
  data$id <- if (is.null(x = group.by)) {
    Idents(object = seurat_obj)[cells, drop = TRUE]
  } else {
    seurat_obj[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data$id)) {
    data$id <- factor(x = data$id)
  }
  data$id <- as.vector(x = data$id)
  
  data$Cell <- rownames(data)
  data <- reshape2::melt(data, variable.name  = "Feature")
  
  #preserve identity order
  if (group.by == "cell.type.ident.by.data.set"){
    data$id <- factor(data$id, levels = levels(seurat_obj$cell.type.ident.by.data.set))
  }else if (group.by == "data.set"){
    data$id <- factor(data$id, levels = levels(seurat_obj$data.set))
  }else if (group.by == "seurat_clusters"){
    data$id <- factor(data$id, levels = levels(seurat_obj$seurat_clusters))
  }else{
    data$id <- factor(data$id, levels = levels(seurat_obj$cell.type.ident))
  }
  
  g <- ggplot(data, aes(Cell, Feature,fill= value)) +
    geom_tile(height = .95, width = 2) +
    scale_fill_distiller(
      palette = "RdYlBu") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y.right = element_text(size=13),
          panel.spacing = unit(.25, "lines"),
          strip.text.x  = element_text(angle = 90, 
                                       vjust = 0.5, 
                                       hjust=.5,
                                       size = 8)) + 
    facet_grid( ~ id, space = 'free', scales = 'free') +
    ylim(rev(levels(data$Feature))) 
  
  return(g)
  
}

# ==== example
# hmap(seurat_obj = seurat_obj,genes = "atoh1a")

volcano_plot <- function(marker_tbl,
                         logFC_hline = c(-1.0, 1.0),
                         lower_bound_logFC = -1.0,
                         upper_bound_logFC = 1.0,
                         mycolors = c("blue", "red", "black"),
                         ident = NULL,
                         gene_col_name = "gene",
                         label.all = FALSE){
  # add a column of NAs
  marker_tbl$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  marker_tbl$diffexpressed[marker_tbl$avg_logFC > upper_bound_logFC & marker_tbl$p_val < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  marker_tbl$diffexpressed[marker_tbl$avg_logFC < lower_bound_logFC & marker_tbl$p_val < 0.05] <- "DOWN"
  
  # 1. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
  names(mycolors) <- c("DOWN", "UP", "NO")
  
  # Now write down the name of genes beside the points...
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  marker_tbl$label <- NA
  labels <- marker_tbl[,gene_col_name][marker_tbl$diffexpressed != "NO"]
  marker_tbl$label[marker_tbl$diffexpressed != "NO"] <- as.character(labels)
  
  if(label.all){
    marker_tbl$label <- as.character(marker_tbl[,gene_col_name])
  }
  
  # plot adding up all layers we have seen so far
  vp <- ggplot(data=marker_tbl, aes(x=avg_logFC, y=-log10(p_val), 
                                    col=diffexpressed, label=label)) +
    geom_point() + 
    theme_minimal() +
    ggrepel::geom_text_repel() +
    scale_color_manual(values=mycolors) +
    geom_vline(xintercept=logFC_hline, col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  if(is.null(ident) == TRUE){
    vp <- vp + labs(title = "Volcano Plot", 
                    subtitle = paste0("Analysis: ", 
                                      ident))
  }else{
    vp <- vp + labs(title = "Volcano Plot", 
                    subtitle = paste0("Analysis: ", 
                                      ident))
  }
  vp <- vp + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(vp)
}

# ==== gse barplot ====
#barplot function is not available for BiocManager V3.10 but appears in latest model
gse_barplot <- function(enrich_obj, showCategory = 10, x_axis = "GeneRatio"){
  #specify whether to plot by counts of genes or ratio of gene counts/setSize
  if (x_axis == "geneRatio" || x_axis == "GeneRatio") {
    x_axis <- "GeneRatio"
  }
  else if (x_axis == "count" || x_axis == "Count") {
    x_axis <- "count"
  }
  
  plot_df <- enrich_obj@result
  
  #showCategory <- 40
  
  #calculate gene counts for core_enrichment 
  gene_count<- plot_df %>% 
    group_by(ID) %>% 
    summarise(count = sum(stringr::str_count(core_enrichment, "/")) + 1) %>%
    mutate(count = as.integer(count))
  
  plot_df <- left_join(plot_df, gene_count, 
                       by = c("ID")) %>% 
    mutate(GeneRatio = as.integer(count)/setSize)
  
  plot_df <- plot_df[order(plot_df$GeneRatio,
                           plot_df$count,decreasing = TRUE),]
  
  #label activated/suppressed terms
  plot_df$group <- "Activated"
  if(any(plot_df$NES <0)){
    #NES measures up and down reg genes
    plot_df[plot_df$NES < 0, ]$group <- "Suppressed"  
  }
  
  #subset top showCategory
  plot_df_sub <- plot_df %>% group_by(group) %>% dplyr::slice(1:showCategory)
  
  plot_df_sub$Description <- factor(plot_df_sub$Description,
                                    levels = plot_df_sub$Description)
  
  #aes_string allows aes to read in character values passed from parameters
  p <-ggplot(data = plot_df_sub,
             aes_string(x = x_axis, y = "Description")) +
    geom_bar(stat='identity', aes(fill = p.adjust)) +
    scale_fill_gradient(low = "purple", high = "red") + 
    facet_grid( ~ group) +
    scale_y_discrete(limits = rev(levels(plot_df_sub$Description)))
  
  
  return(p)
  
}


# ==== plot all enrichment graphs from clusterProfiler ====
#plot all standard plot output from GSEA analysis
make_enrich_plot <- function(enrich_obj,
                             marker_tbl, 
                             dir_name,
                             dir_name_sub, 
                             showCategory = 10, 
                             analysis.type = c("GSEA", "KEGG", "REACTOME"),
                             is.enrich.reactome = F , 
                             save = T,
                              color_by = "p.adjust",
                              x_axis = "count") {
  
  require(DOSE)
  
  if(analysis.type == "GSEA"){
    title <- as.character("Gene Set Enrichment")
  }else if(analysis.type == "KEGG" | analysis.type == "REACTOME"){
    title <- as.character("Enriched Pathways")
  }
  
  if(is.enrich.reactome){
    rh <- enrichplot::heatplot(enrich_obj$reactome_obj,
                   foldChange = enrich_obj$reactome_gene_list,
                   showCategory = showCategory)
    rh <- rh + scale_fill_distiller(
      palette = "RdYlBu")

    #append our gene symbol to x axis label, replace entrez ID
    rh_data <- left_join(rh$data,
                         marker_tbl[,c("Gene.name.uniq","ENTREZID")],
                         by = c("Gene" = "ENTREZID"))

    rh <- rh + scale_x_discrete(labels = factor(rh_data$Gene.name.uniq)) +
       labs(title = title,
             subtitle = paste0("Analysis: ",
                               unique(marker_tbl$ident)),
             x = "enrichment distribution")
    
    d <- NULL
    
    r <- NULL
    
    if(save){
      png(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_heatplot_",dir_name_sub, ".png")),
          width = 14,
          height = 8,
          units = "in",
          res = 300)
      print(rh)
      dev.off()
    }
    
  }else{
    rh <- NULL
    
    d <- enrichplot::dotplot(enrich_obj, 
                             showCategory=showCategory, 
                             split=".sign", color = color_by) + 
      facet_grid(.~.sign)
    
    d <- d + labs(title = title, 
                  subtitle = paste0("Analysis: ", 
                                    unique(marker_tbl$ident)))
    
    r <- enrichplot::ridgeplot(enrich_obj, showCategory =showCategory,
                               fill = color_by) 
    r <- r + labs(title = title, 
                  subtitle = paste0("Analysis: ", 
                                    unique(marker_tbl$ident)),
                  x = "enrichment distribution")
    
    b <- gse_barplot(enrich_obj = enrich_obj,
                     showCategory=showCategory,
                     x_axis = x_axis)

    b <- b + labs(title = title,
                 subtitle = paste0("Analysis: ",
                                   unique(marker_tbl$ident)))
    
    if(save){
      png(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_dotplot_",dir_name_sub, ".png")),
          width = 14,
          height = 8,
          units = "in",
          res = 300)
      print(d)
      dev.off()
  
      png(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_ridgeplot_",dir_name_sub, ".png")),
          width = 14,
          height = 8,
          units = "in",
          res = 300)
      print(r)
      dev.off()
      
      png(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_barplot_",dir_name_sub, ".png")),
          width = 14,
          height = 8,
          units = "in",
          res = 300)
      print(b)
      dev.off()
      
      
    }
    
  }

  return(list(dotplot = d, ridgeplot = r, barplot = b, heatplot = rh))

}

# make expression line plots by timepoints

line_plots_by_time <- function(seurat_obj, genes, group.by = "samples"){
  require(reshape2)
  
  plotting_df <- 
    SeuratExtensions::get_expression_tbl(seurat_obj = seurat_obj,
                                         genes = genes,
                                         idents = group.by, 
                                         slot = "data"
    )
  
  plotting_df <- melt(plotting_df)
  
  names(plotting_df)[names(plotting_df) == "variable"] <- "ident"
  names(plotting_df)[names(plotting_df) == "value"] <- "expression"
  
  plotting_df$data.set <- gsub("(.+?)(\\_.*)", "\\1",plotting_df$ident)
  plotting_df$data.set <- factor(plotting_df$data.set, 
                                 levels = levels(obj_integrated$data.set))
  
  plotting_df$Treatment <- gsub(".*_", "\\1",plotting_df$ident)
  plotting_df$Treatment <- factor(plotting_df$Treatment, 
                                  levels = levels(seurat_obj$treatment))
  
  
  plotting_df$scaled.exp <-  scale(plotting_df$expression)
  
  lp <- ggplot(data=plotting_df,
               aes(x=data.set,
                   y=scaled.exp, 
                   colour=Treatment,  
                   group = Treatment)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ Gene.name.uniq, ncol = 2) +
    ylab(label = "Scaled Average Expression") + 
    xlab(label = "Time") +
    theme(strip.text = element_text(face = "italic"))
  
  return(lp)
  
}
