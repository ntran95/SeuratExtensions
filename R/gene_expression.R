get_expression_tbl <- function(seurat_obj, genes, idents, slot = "data"){
  Idents(seurat_obj) <- idents
  
  avg.mtx <-  Seurat::AverageExpression(object = seurat_obj,
                                        features = genes,
                                        assays = "RNA",
                                        slot = slot)$RNA
  
  if(slot == "scale.data"){
    avg.mtx <- Seurat::MinMax(avg.mtx,min = -2.5, max = 2.5)
  }
  
  if(!exists("gene_info")){
    gene_info <- read.delim(paste0("/",volumes,"/projects/nt2473/Analysis/Data/",
                                   "gene-lists/Danio_Features_unique_Ens91_v2.tsv"),
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  
  avg.mtx$Gene.name.uniq <- rownames(avg.mtx)
  
  expression_tbl <- dplyr::right_join(gene_info,avg.mtx,by = "Gene.name.uniq")
  
  #match order with genes vector
  expression_tbl<- expression_tbl[match(genes, expression_tbl$Gene.name.uniq),]
  
  return(expression_tbl)
}