corrleatedGenes <- function(seurat_obj, goi,
  cell_type = unique(Idents(seurat_obj)),
  condition = "homeo", n_results = 200, verbose = TRUE) {
  
  if (length(goi) > 1) {
    stop("Please enter a single gene")
  }

  mat <- seurat_obj@assays$RNA@data
  mat_raw <- as.matrix(mat)
  
  mat_sub <- mat_raw[, seurat_obj[["data.set"]][,1] %in% condition &
    seurat_obj[['cell.type.ident']][,1] %in% cell_type]

  gene <- as.numeric(mat_sub[goi,])
  cor_vals <- apply(mat_sub, 1, function(x) {cor(gene,x)})
  cor_vals <- sort(cor_vals, decreasing = TRUE)
  cor_vals <- c(head(cor_vals, n_results), tail(cor_vals, n_results))
  
  cor_df <- data.frame(Gene.name.uniq = names(cor_vals),
      correlation_coeff = cor_vals)

  if (verbose) {cat(rownames(cor_df),sep="\n")}
  return(cor_df)
}
