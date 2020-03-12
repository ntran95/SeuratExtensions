evaluate_cell_selection <- FALSE
if (evaluate_cell_selection) {

# Other normalization methods "RC" and "CLR" "LogNormalize"
obj_test <- obj_integrated
DefaultAssay(obj_test) <- "RNA"
obj_test <- NormalizeData(obj_test,
  normalization.method = "CLR", verbose = TRUE)

# ----- Get BCs for cell type
# order by cell type within treatment
meta <- obj_test@meta.data
meta$cell.type.and.trt <- paste0(meta$cell.type.ident, "_", meta$data.set)

to_order <- unique(meta$cell.type.and.trt)
cell_nm <- sort(unique(meta$cell.type.ident))

# ! Change for project or data set
cell_nm <- cell_nm[c(1,3,5,2,4,6:13)]
trt <- sort(unique(meta$data.set))

# Remove unwanted cell types
rm_cell_grps <- FALSE
if (rm_cell_grps) {
  cell_nm <- cell_nm[-c(9:14)]
  meta <- meta[meta$cell.type.ident %in% cell_nm,]
}

dim(meta)
ordered_list <- list()

# Order each cell type within a timpeoint, or timepoint wihtin cell type
cell_type_within_time <- FALSE
if (cell_type_within_time) {
  for (i in seq_along(trt)) {
    to_match <- sort(to_order[grep(paste0("_", trt[i]), to_order)])
    ordered_list[[i]] <- to_match[match(cell_nm, gsub("_.*","", to_match))]
  }
} else {
  to_order <- sort(to_order)
  for (i in seq_along(cell_nm)) {
    to_match <- sort(to_order[grep(cell_nm[i], to_order)])
    ordered_list[[i]] <- to_match[match(trt, gsub(".*_","", to_match))]
  }
}

# Collapse ordered list
ordered_vec <- unlist(ordered_list); ordered_vec

# ordering meta data to new cell type arrangement
meta$cell.type.and.trt <- factor(
  meta$cell.type.and.trt, levels = ordered_vec, ordered = TRUE)
meta <- meta[order(meta$cell.type.and.trt),]

# BCs are now in order, subset them
ordered_BCs <- rownames(meta)


# ====== Subset expression matrix
scaled <- TRUE
if (scaled) {
  DefaultAssay(obj_test) <- "RNA"
  obj_test <- ScaleData(obj_test)
  exp_mat <- obj_test@assays$RNA@scale.data
} else {
  exp_mat <- obj_test@assays$RNA@data
}
exp_mat <- as.matrix(exp_mat)

mat_sub <- exp_mat[,colnames(exp_mat) %in% rownames(meta)]
mat_sub <- mat_sub[,match(rownames(meta), colnames(mat_sub))]

# Check that meta data BCs are in order with expression matrix
head(colnames(mat_sub)); head(rownames(meta))

# ====== Avg cell types and timepoints
# Get factor level index changes for treatment
cell_type_idx <- match(unique(meta$cell.type.and.trt), meta$cell.type.and.trt)
cell_type_idx

# Create list of dfs with average exp
trt_avg_list <- lapply(seq_along(cell_type_idx), function(i) {
  if (i < length(cell_type_idx)) {
    trt_ind <- (cell_type_idx[i]) : (cell_type_idx[i + 1] - 1)
  } else {
    trt_ind <- cell_type_idx[i] : ncol(mat_sub)}

  trt_avg_exp <- apply(mat_sub[,trt_ind], 1, mean)
  return(data.frame(trt_avg_exp))
})

length(trt_avg_list)
trt_avg_df <- bind_cols(trt_avg_list)
dim(trt_avg_df)

# Get fold change
if (scaled) {
  logFC_df0 <- trt_avg_df
  } else {
  logFC_df0 <- as.data.frame(sapply(c(1:ncol(trt_avg_df)), function(i) {
      (trt_avg_df[,i]) - (trt_avg_df[,1])
  }))
}

dim(logFC_df0)
head(logFC_df0[2:10,1])

rownames(logFC_df0) <- rownames(mat_sub)
colnames(logFC_df0) <- unique(meta$cell.type.and.trt)

} # End Select Cells

# ====== Get genes of interest
cell_specific <- TRUE
if (cell_specific) {
  file_name <- "MarkerList_specific_"
  all_markers <- readRDS(dataPath(paste0(file_name, script_name,"_.RDS")))
  folder <- "diff_exp_cell_specific/"
} else {
  file_name <- "MarkerList_"
  all_markers <- readRDS(dataPath(paste0(file_name, script_name,"_.RDS")))
  folder <- "diff_exp_by_cell_type/"
}

goi <- TRUE
if (goi) {
  unique(all_markers$cell.type.and.trt)
  marker_sub <- all_markers[
    all_markers$cell.type.and.trt %in% c(
      "central-cells_1hr", "central-cells_3hr", "central-cells_5hr"),]
  to_plot <- marker_sub[order(marker_sub$avg_logFC), "Gene.name.uniq"]
  to_plot <- marker_sub$Gene.name.uniq[1:length(to_plot)]
  length(to_plot)
} else {
  to_plot <- to_plot[1:length(to_plot)]}; length(to_plot)

# Systemic genes from pbulk all cell, all timepoints
sys_genes <- FALSE
if (sys_genes) {
  diff_table <- readRDS(dataPath(paste0(
  "all_cells_pbulk_by_time_", script_name,"_.RDS")))

  min0 <- diff_table[diff_table$timepoint == "0min", "Gene.name.uniq"][1:25]
  min30 <- diff_table[diff_table$timepoint == "30min", "Gene.name.uniq"][1:25]
  hour1 <- diff_table[diff_table$timepoint == "1hr", "Gene.name.uniq"][1:25]
  sys_genes <- c(min0, min30, hour1)
}

remove_duplicates <- TRUE
remove_genes <- TRUE
if (remove_genes) {
  if (remove_duplicates) {
    to_plot <- to_plot[!duplicated(to_plot)]
  }
  remove_these <- c("aqp3a","si:ch211-152f23.5","prr15la","gadd45ga",
    "sparc","krt18","krt8","krt92","rgcc","id2a", "crabp2a", "mfng")
  to_plot <- to_plot[!(to_plot %in% remove_these)]
  length(to_plot)
}

add_genes <- FALSE
if (add_genes) {
  add_these <- c(sys_genes)
  to_plot <- c(sys_genes, to_plot)
}
length(to_plot)

# ----- Create color annotation for pheatmap
group_names <- colnames(logFC_df0)
trt <- gsub(".*_","", group_names)

anno_column <- data.frame(timepoint = factor(trt,
  levels = levels(meta$data.set), ordered = TRUE))
rownames(anno_column) <- group_names

names(trt_colors) <-  unique(trt)
dataset <- trt_colors

# color variable "timepoint" must be the SAME as the annotation
anno_cols <- list(timepoint = dataset); anno_cols

gaps <- match(unique(trt), trt)
gaps <- (gaps - 1)[-1]

# ----- Slice matrix on genes of interest and plotlistt
goi_mat <- logFC_df0[rownames(logFC_df0) %in% to_plot,]
dim(goi_mat)
goi_mat <- goi_mat[match(to_plot, rownames(goi_mat)),]

change_color <- FALSE
if (change_color) {
library(RColorBrewer)
color_gradient <- colorRampPalette(
  rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
}

# height = 40 for 100 genes
png(figurePath("phmap_test.png"),
  units = "in", res = 200, width = 20, height = 20)
pheatmap::pheatmap(goi_mat, cluster_rows = TRUE, cluster_cols = FALSE,
  color = viridis::viridis(100), annotation_col = NULL, legend = FALSE,
  annotation_colors = anno_cols, gaps_col = c(7,14), annotation_names_col = FALSE,
  annotation_legend = FALSE)
dev.off()
