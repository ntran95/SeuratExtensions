cleanUMAP <- function(plot_obj, dark = FALSE, axis_title_size = 16,
  plot_title = "") {
  if (dark == TRUE) {
    plot_obj <- plot_obj + labs(title = element_text(plot_title), x = "UMAP 1", y = "UMAP 2") +
      theme(plot.background = element_rect(fill = "black",
          colour = "black", size = 0.5, linetype = "solid"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.line.x = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      axis.title = element_text(size = axis_title_size, color = "white"),
      legend.text = element_text(color = c("white"), size = 10)) +
      ggtitle(plot_title)
    } else {
    plot_obj <- plot_obj + labs(title = element_text(plot_title), x = "UMAP 1", y = "UMAP 2") +
      theme(plot.title = element_text(plot_title), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title = element_text(size = axis_title_size))
  }
}
