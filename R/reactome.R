#get_enrich_reactome_obj uses unsorted entrez ID vector, cannot distinguish activated and suppressed terms
#use this to build heatplot()
make_enrich_reactome_obj <- function(marker_tbl,
                             organism = "org.Dr.eg.db",
                             reactome_organism  = "zebrafish",
                             minGSSize = 3, 
                             maxGSSize = 800,
                             pvalueCutoff = 0.05){
  
  #convert ensembl --> entrez (KEGG only takes in entrez)
  entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                      fromType = "ENSEMBL", 
                                      toType = "ENTREZID", 
                                      OrgDb=organism) 
  #remove duplicate IDs
  dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  marker_tbl_entrez <- inner_join(marker_tbl, dedup_ids, 
                                  by = c("Gene.stable.ID" = "ENSEMBL"))
  
  # Create a vector of the gene unuiverse
  reactome_gene_list <- marker_tbl_entrez$avg_logFC
  
  # Name vector with ENTREZ ids
  names(reactome_gene_list) <- marker_tbl_entrez$ENTREZID
  
  # omit any NA values 
  reactome_gene_list<-na.omit(reactome_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  reactome_gene_list <- sort(reactome_gene_list, decreasing = TRUE)
  
  reactome_obj <- ReactomePA::enrichPathway(marker_tbl_entrez$ENTREZID,
                                            organism = reactome_organism,
                                            pvalueCutoff = pvalueCutoff,
                                            pAdjustMethod = "BH",
                                            qvalueCutoff = 0.2,
                                            universe,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize)
  
  return(list(reactome_obj = reactome_obj,
              reactome_gene_list = reactome_gene_list))
}

#get_reactome_obj uses sorted geneList, distinguish between activated and suppressed terms
make_reactome_obj <- function(marker_tbl,
                             organism = "org.Dr.eg.db",
                             reactome_organism  = "zebrafish",
                             minGSSize = 3, 
                             maxGSSize = 800,
                             pvalueCutoff = .05){
  #convert ensembl --> entrez (KEGG only takes in entrez)
  entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                      fromType = "ENSEMBL", 
                                      toType = "ENTREZID", 
                                      OrgDb=organism) 
  #remove duplicate IDs
  dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  marker_tbl_entrez <- inner_join(marker_tbl, dedup_ids, 
                                  by = c("Gene.stable.ID" = "ENSEMBL"))
  
  # Create a vector of the gene unuiverse
  reactome_gene_list <- marker_tbl_entrez$avg_logFC
  
  # Name vector with ENTREZ ids
  names(reactome_gene_list) <- marker_tbl_entrez$ENTREZID
  
  # omit any NA values 
  reactome_gene_list<-na.omit(reactome_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  reactome_gene_list <- sort(reactome_gene_list, decreasing = TRUE)
  
  set.seed(123)
  reactome <- ReactomePA::gsePathway(geneList = reactome_gene_list,
                                     organism = reactome_organism,
                                     nPerm = 10000,
                                     minGSSize = minGSSize,
                                     maxGSSize = maxGSSize,
                                     seed = T,
                                     pvalueCutoff = pvalueCutoff)
  
  return(reactome)
}

add_reactome_to_marker_tbl <- function(marker_tbl, 
                                   reactome,
                                   organism = "org.Dr.eg.db",
                                   core_enrichment = TRUE){
  #convert ensembl --> entrez (reactome only takes in entrez)
  entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                      fromType = "ENSEMBL", 
                                      toType = "ENTREZID", 
                                      OrgDb=organism) 
  #remove duplicate IDs
  dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  marker_tbl_entrez <- left_join(marker_tbl, dedup_ids, 
                                 by = c("Gene.stable.ID" = "ENSEMBL"))
  
  reactome_results <- reactome@result
  
  if(core_enrichment){
  #expand ensembl ID to individual rows in "core_enrichment" column
  reactome_results_expand <- reactome_results %>% 
    mutate(core_enrichment = strsplit(as.character(core_enrichment), "/")) %>%
    tidyr::unnest(core_enrichment)
  
  #subset for relevant info
  reactome_results_expand <- reactome_results_expand[,c("ID",
                                                "Description", 
                                                "core_enrichment")]
  
  #collapse go terms by ensembl ID
  reactome_results_collapse <- reactome_results_expand %>% 
    group_by_at(vars(core_enrichment)) %>%
    summarise_at(vars(ID,Description), paste, collapse = ",")
  
  #rename columns
  reactome_cols <- 
    paste0("reactome.", colnames(reactome_results_collapse[,c(2:3)]))
  
  colnames(reactome_results_collapse)[2:3] <-reactome_cols
  
  
  marker_tbl_reactome <- left_join(marker_tbl_entrez, 
                               reactome_results_collapse,
                               by= c( "ENTREZID" = "core_enrichment"))
  }else{ #if using enrichPathway()
    #expand ensembl ID to individual rows in "geneID" column
    reactome_results_expand <- reactome_results %>% 
      mutate(geneID = strsplit(as.character(geneID), "/")) %>%
      tidyr::unnest(geneID)
    
    #subset for relevant info
    reactome_results_expand <- reactome_results_expand[,c("ID",
                                                          "Description", 
                                                          "geneID")]
    
    #collapse go terms by ensembl ID
    reactome_results_collapse <- reactome_results_expand %>% 
      group_by_at(vars(geneID)) %>%
      summarise_at(vars(ID,Description), paste, collapse = ",")
    
    #rename columns
    reactome_cols <- 
      paste0("reactome.", colnames(reactome_results_collapse[,c(2:3)]))
    
    colnames(reactome_results_collapse)[2:3] <-reactome_cols
    
    
    marker_tbl_reactome <- left_join(marker_tbl_entrez, 
                                     reactome_results_collapse,
                                     by= c( "ENTREZID" = "geneID"))
  }
  
  return(marker_tbl_reactome)
  
}
