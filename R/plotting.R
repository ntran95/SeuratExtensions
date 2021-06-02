hmap <- function(seurat_obj, genes, ident = "cell.type.ident"){
  dotplot <- Seurat::DotPlot(seurat_obj, features = genes,
                     group.by = ident)
  
  g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, 
                                width = 2, height = 2)) +
    geom_tile(color = "gray", size = 1) +
    scale_fill_distiller(
      palette = "RdYlBu") +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5, 
                                     hjust=.5,
                                     size = 13),
          axis.title.y.right = element_text(size=13),
          strip.text.x  = element_text(vjust = 0.5, 
                                       hjust=.5,
                                       size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
  scale_x_discrete(expand=c(0,0)) + #expand tile 
    scale_y_discrete(expand=c(0,0))+ #expand tile +
  ylim(rev(levels(dotplot$data$features.plot))) 
  
  if(ident == "cell.type.ident.by.data.set"){
    dotplot$data$groupIdent <- gsub("(.+?)(\\_.*)", "\\1",dotplot$data$id)
    dotplot$data$groupIdent <- factor(dotplot$data$groupIdent,levels=levels(seurat_obj$cell.type.ident))
    
    g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) +
      geom_tile(color = "gray", size = 1) +
      scale_fill_distiller(
        palette = "RdYlBu") +
      hrbrthemes::theme_ipsum(base_family = "sans") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
            axis.title.y.right = element_text(size=13),panel.spacing = unit(.35, "lines"),
            strip.text.x  = element_text(vjust = 0.5, hjust=.5,size = 12)) +
      facet_grid( ~ groupIdent, scales='free_x')+
      ylim(rev(levels(dotplot$data$features.plot))) +
      theme_ipsum(base_family = "sans")

  }
  
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
  vp <- ggplot(data=marker_tbl, aes(x=avg_logFC, y=-log10(p_val), col=diffexpressed, label=label)) +
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

#plot all standard plot output from GSEA analysis
make_enrich_plot <- function(enrich_obj,
                             marker_tbl, 
                             dir_name,
                             dir_name_sub, 
                             showCategory = 10, 
                             analysis.type = c("GSEA", "KEGG", "REACTOME"),
                             is.enrich.reactome = F , 
                             save = T,
                            color_by = "p.adjust") {
  
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
      pdf(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_heatplot_",dir_name_sub, ".pdf")),
          width = 14,
          height = 8)
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
    
    if(save){
      pdf(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_dotplot_",dir_name_sub, ".pdf")),
          width = 14,
          height = 8)
      print(d)
      dev.off()
  
      pdf(figurePath(paste0(dir_name, "/", dir_name_sub, "/",
                            analysis.type,"_ridgeplot_",dir_name_sub, ".pdf")),
          width = 14,
          height = 8)
      print(r)
      dev.off()
    }
    
  }

  return(list(dotplot = d, ridgeplot = r, heatplot = rh))

}
