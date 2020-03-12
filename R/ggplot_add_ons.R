
cleanUMAP <- function(dark = FALSE){
  if (dark == TRUE) {
    p1 <- p1 + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), plot.background = element_rect(
        fill = "black", colour = "black", size = 0.5, linetype = "solid"),
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), axis.title = element_text(size = 16,
        color = "white"), legend.text = element_text(color = c("white"),
        size = 10))
    } else {
    p1 <- p1 + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),  axis.title = element_text(size = 16))
  }
}