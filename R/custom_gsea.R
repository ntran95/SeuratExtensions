# read in user annotated custom terms for universal enrichment
# pools terms from multiple databases into single database

#use GSEA() to get activated and suppressed terms in dotplot, 
#reads in ranked geneList
make_custom_GSEA <- function(marker_tbl,
                             organism = "org.Dr.eg.db", #zf annotations
                             minGSSize = 2, 
                             maxGSSize = 800,
                             pvalueCutoff = 0.05,
                             addEntrezID = F, #does the ENTREZID column need to be appended to marker_tbl?
                             TERM2GENE=term2gene, 
                             TERM2NAME=term2name){
  
  if(addEntrezID){
  #convert ensembl --> entrez (custom db is in entrez ID due to kegg end reactome)
  entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                      fromType = "ENSEMBL", 
                                      toType = "ENTREZID", 
                                      OrgDb=organism) 
  
  #remove duplicate IDs
  dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  marker_tbl <- inner_join(marker_tbl, dedup_ids, 
                                  by = c("Gene.stable.ID" = "ENSEMBL"))
  }
  
  #filter out rows with NA in entrez ID
  marker_tbl <- marker_tbl[!is.na(marker_tbl$ENTREZID),]
  
  # Create a vector of the gene 
  gene_list <- marker_tbl$avg_logFC
  
  # Name vector with ENTREZ ids
  names(gene_list) <- marker_tbl$ENTREZID
  
  # omit any NA values 
  gene_list<-na.omit(gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  set.seed(123)
  
  custom_gsea <- clusterProfiler::GSEA(gene_list, 
                        TERM2GENE=term2gene, 
                        TERM2NAME=term2name,
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        pvalueCutoff = pvalueCutoff,
                        nPerm = 10000,
                        seed = T)
  
  return(custom_gsea)
}

make_custom_enricher <- function(marker_tbl,
                                  organism = "org.Dr.eg.db", #zf annotations
                                  minGSSize = 2, 
                                  maxGSSize = 800,
                                  pvalueCutoff = 0.05,
                                  addEntrezID = F, #does the ENTREZID column need to be appended to marker_tbl?
                                  TERM2GENE=term2gene, 
                                  TERM2NAME=term2name){
  
  if(addEntrezID){
    #convert ensembl --> entrez (custom db is in entrez ID due to kegg end reactome)
    entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                        fromType = "ENSEMBL", 
                                        toType = "ENTREZID", 
                                        OrgDb=organism) 
    
    #remove duplicate IDs
    dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
    
    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    marker_tbl <- inner_join(marker_tbl, dedup_ids, 
                             by = c("Gene.stable.ID" = "ENSEMBL"))
  }
  
  #filter out rows with NA in entrez ID
  marker_tbl <- marker_tbl[!is.na(marker_tbl$ENTREZID),]
  
  gene_list <- marker_tbl$ENTREZID
  
  custom_gsea <- clusterProfiler::enricher(gene_list, 
                                       TERM2GENE=term2gene, 
                                       TERM2NAME=term2name,
                                       minGSSize = minGSSize,
                                       maxGSSize = maxGSSize,
                                       pvalueCutoff = pvalueCutoff)
  
  return(custom_gsea)
}

add_custom_gsea_to_marker_tbl <- function(marker_tbl, 
                                          custom_gsea,
                                          addEntrezID = F,
                                          combined_term_df) {#does the ENTREZID column need to be appended to marker_tbl?){
   if(addEntrezID){
  #convert ensembl --> entrez (custom db is in entrez ID due to kegg end reactome)
  entrez_ids <- clusterProfiler::bitr(marker_tbl$Gene.stable.ID, 
                                      fromType = "ENSEMBL", 
                                      toType = "ENTREZID", 
                                      OrgDb=organism) 
  
  #remove duplicate IDs
  dedup_ids = entrez_ids[!duplicated(entrez_ids[c("ENSEMBL")]),] 
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  marker_tbl <- inner_join(marker_tbl, dedup_ids, 
                                  by = c("Gene.stable.ID" = "ENSEMBL"))
   }
  
  gsea_results <- custom_gsea@result
  
  #expand ensembl ID to individual rows in "core_enrichment" column
  results_expand <- gsea_results %>% 
    mutate(core_enrichment = strsplit(as.character(core_enrichment), "/")) %>%
    tidyr::unnest(core_enrichment)
  
  #subset for relevant info
  results_expand <- results_expand[,c("ID",
                                      "Description", 
                                      "core_enrichment")]
  
  results_expand <- inner_join(results_expand,
                              combined_term_df[,c("ID","database")],
                              by = "ID")
  
  #collapse go terms by ensembl ID
  results_collapse <- results_expand %>% 
    group_by_at(vars(core_enrichment)) %>%
    summarise_at(vars(ID,Description,database), paste, collapse = ",")
  
  # #rename columns
  # kegg_cols <- 
  #   paste0("KEGG.", colnames(kegg_results_collapse[,c(2:3)]))
  # 
  # colnames(kegg_results_collapse)[2:3] <-kegg_cols
  # 
  
  marker_tbl <- left_join(marker_tbl, 
                          results_collapse,
                          by= c( "ENTREZID" = "core_enrichment"))
  
  return(marker_tbl)
}
  