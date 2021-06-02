make_kegg_obj <- function(marker_tbl,
                          organism = "org.Dr.eg.db",
                          kegg_organism  = "dre",
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
  kegg_gene_list <- marker_tbl_entrez$avg_logFC
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- marker_tbl_entrez$ENTREZID
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
  
  #create gseKEGG object
  #kegg_organism <- "dre"
  set.seed(123)
  kegg <- clusterProfiler::gseKEGG(geneList = kegg_gene_list,
                  organism     = kegg_organism,
                  nPerm        = 10000,
                  minGSSize    = minGSSize,
                  maxGSSize    = maxGSSize,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = "none",
                  keyType       = "ncbi-geneid",
                  seed = T)
  
  return(kegg)
}


#get_enrich_kegg_obj uses unsorted entrez ID vector, cannot distinguish activated and suppressed terms
#use this to build heatplot()
make_enrich_kegg_obj <- function( marker_tbl,
                                  organism = "org.Dr.eg.db",
                                  kegg_organism  = "dre",
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
  #kegg_gene_list <- marker_tbl_entrez$avg_logFC
  
  # Name vector with ENTREZ ids
  #names(kegg_gene_list) <- marker_tbl_entrez$ENTREZID
  
  # omit any NA values 
 # kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  #kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
  
  kegg_obj <- clusterProfiler::enrichKEGG(marker_tbl_entrez$ENTREZID,
                                            organism = kegg_organism,
                                            pvalueCutoff = pvalueCutoff,
                                            pAdjustMethod = "BH",
                                            qvalueCutoff = 0.2,
                                            universe,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize,
                                            keyType= "ncbi-geneid")
  
  # return(list(kegg_obj = kegg_obj,
  #             kegg_gene_list = kegg_gene_list))
  # 
  return(kegg_obj = kegg_obj)
}


add_kegg_to_marker_tbl <- function(marker_tbl, 
                                   kegg,
                                   organism = "org.Dr.eg.db",
                                   kegg_organism  = "dre",
                                   core_enrichment = TRUE){
  #convert ensembl --> entrez (KEGG only takes in entrez)
  entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                      fromType = "ENSEMBL", 
                                      toType = "ENTREZID", 
                                      OrgDb=organism) 
  #remove duplicate IDs
  dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  marker_tbl_entrez <- left_join(marker_tbl, dedup_ids, 
                                  by = c("Gene.stable.ID" = "ENSEMBL"))
  
  kegg_results <- kegg@result
  
  if(core_enrichment){
    #expand ensembl ID to individual rows in "core_enrichment" column
    kegg_results_expand <- kegg_results %>% 
      mutate(core_enrichment = strsplit(as.character(core_enrichment), "/")) %>%
      tidyr::unnest(core_enrichment)
    
    #subset for relevant info
    kegg_results_expand <- kegg_results_expand[,c("ID",
                                                  "Description", 
                                                  "core_enrichment")]
    
    #collapse go terms by ensembl ID
    kegg_results_collapse <- kegg_results_expand %>% 
      group_by_at(vars(core_enrichment)) %>%
      summarise_at(vars(ID,Description), paste, collapse = ",")
    
    #rename columns
    kegg_cols <- 
      paste0("KEGG.", colnames(kegg_results_collapse[,c(2:3)]))
    
    colnames(kegg_results_collapse)[2:3] <-kegg_cols
    
    
    marker_tbl_kegg <- left_join(marker_tbl_entrez, 
                                 kegg_results_collapse,
                                 by= c( "ENTREZID" = "core_enrichment"))
  }else{ #if using enrichKEGG()
    #expand ensembl ID to individual rows in "geneID" column
    kegg_results_expand <- kegg_results %>% 
      mutate(geneID = strsplit(as.character(geneID), "/")) %>%
      tidyr::unnest(geneID)
    
    #subset for relevant info
    kegg_results_expand <- kegg_results_expand[,c("ID",
                                                  "Description", 
                                                  "geneID")]
    
    #collapse go terms by ensembl ID
    kegg_results_collapse <- kegg_results_expand %>% 
      group_by_at(vars(geneID)) %>%
      summarise_at(vars(ID,Description), paste, collapse = ",")
    
    #rename columns
    kegg_cols <- 
      paste0("KEGG.", colnames(kegg_results_collapse[,c(2:3)]))
    
    colnames(kegg_results_collapse)[2:3] <-kegg_cols
    
    
    marker_tbl_kegg <- left_join(marker_tbl_entrez, 
                                 kegg_results_collapse,
                                 by= c( "ENTREZID" = "geneID"))
  }
  return(marker_tbl_kegg)
  
}
