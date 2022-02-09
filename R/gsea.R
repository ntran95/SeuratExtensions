make_gsea_obj <- function(marker_tbl,
                          organism = "org.Dr.eg.db",
                          minGSSize = 3, 
                          maxGSSize = 800,
                          pvalueCutoff = 0.05){
  
  library(organism, character.only = TRUE)

  # we want the log2 fold change 
  original_gene_list <- marker_tbl$avg_logFC

  # name the vector
  names(original_gene_list) <- marker_tbl$Gene.stable.ID

  # omit any NA values 
  gene_list<-na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)

  set.seed(123)
  gse <- clusterProfiler::gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = minGSSize, 
             maxGSSize = maxGSSize, 
             pvalueCutoff = pvalueCutoff, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none",
             seed = T)
  
  return(gse)
}

#get_enrich_gsea_ obj uses unsorted entrez ID vector, cannot distinguish activated and suppressed terms
#use this to build heatplot()
make_enrich_gsea_obj <- function( marker_tbl,
                                  organism = "org.Dr.eg.db",
                                  minGSSize = 3, 
                                  maxGSSize = 800,
                                  keyType= "ENSEMBL",
                                  pvalueCutoff = 0.05){
  
  library(organism, character.only = TRUE)
  
  # we want the log2 fold change 
  original_gene_list <- marker_tbl$avg_logFC
  
  # name the vector
  names(original_gene_list) <- marker_tbl$Gene.stable.ID
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse_obj <- clusterProfiler::enrichGO(gene = names(gene_list),
                                        OrgDb = organism,
                                        keyType = keyType,
                                        ont = "ALL",
                                        pvalueCutoff = pvalueCutoff,
                                        pAdjustMethod = "BH",
                                        universe,
                                        qvalueCutoff = 0.2,
                                        minGSSize = minGSSize,
                                        maxGSSize = maxGSSize)
  
  return(list(gse_obj = gse_obj,
              gene_list = gene_list))
}

add_gsea_to_marker_tbl <- function(marker_tbl, 
                                   gsea, 
                                   core_enrichment =TRUE){
  
  #append Go Terms to marker_tbl
  gse_results <- gsea@result
  if (core_enrichment){
    #expand ensembl ID to individual rows in "core_enrichment" column
    gse_results_expand <- gse_results %>% 
      mutate(core_enrichment = strsplit(as.character(core_enrichment), "/")) %>%
      tidyr::unnest(core_enrichment)
    
    colnames(gse_results_expand)
    #subset for relevant infp
    gse_results_expand <- gse_results_expand[,c("ID",
                                                "Description", 
                                                "core_enrichment")]
    #collapse go terms by ensembl ID
    gse_results_collapse <- gse_results_expand %>% 
      group_by_at(vars(core_enrichment)) %>%
      summarise_at(vars(ID,Description), paste, collapse = ",")
    
    #rename columns
    gsea_cols <- 
      paste0("GSEA.", colnames(gse_results_collapse[,c(2:3)]))
    
    colnames(gse_results_collapse)[2:3] <- gsea_cols
    
    marker_tbl_gse <- left_join(marker_tbl, 
                                gse_results_collapse,
                                by= c( "Gene.stable.ID" = "core_enrichment"))
  }else{ #if using enrichGO()
    #expand ensembl ID to individual rows in "core_enrichment" column
    gse_results_expand <- gse_results %>% 
      mutate(geneID = strsplit(as.character(geneID), "/")) %>%
      tidyr::unnest(geneID)
    
    colnames(gse_results_expand)
    #subset for relevant infp
    gse_results_expand <- gse_results_expand[,c("ID",
                                                "Description", 
                                                "geneID")]
    #collapse go terms by ensembl ID
    gse_results_collapse <- gse_results_expand %>% 
      group_by_at(vars(geneID)) %>%
      summarise_at(vars(ID,Description), paste, collapse = ",")
    
    #rename columns
    gsea_cols <- 
      paste0("GO", colnames(gse_results_collapse[,c(2:3)]))
    
    colnames(gse_results_collapse)[2:3] <- gsea_cols
    
    marker_tbl_gse <- left_join(marker_tbl, 
                                gse_results_collapse,
                                by= c( "Gene.stable.ID" = "geneID"))
  }
  
  return(marker_tbl_gse)

}
