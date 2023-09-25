### Coefficient heatmap generator
heatmap_cells_coefficient <- function(seuratObject, 
                                      store_dir, 
                                      store_folder = "Heatmap_of_the_coefficient",
                                      cell_ids,
                                      plot_by_clusters = TRUE,
                                      figureName = "cells",
                                      reduction_pattern = "inmf") {
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # Let's get the coefficient matrix
  coef = as.data.frame(seuratObject@reductions[[reduction_pattern]]@cell.embeddings)
  colnames(coef) = str_c("GEP_", parse_number(colnames(coef)))
  
  coef = order_factorized_matrix(coef)
  
  cells_order = c()
  
  for (i in levels(seuratObject)){
    cell_idents = WhichCells(seuratObject, idents = i)
    cells_order = c(cells_order, cell_idents)
  }
  
  if (plot_by_clusters){
    mat_H = coef[cells_order, ]
  }
  else {
    mat_H = mat_H
  }
  
  mat_H = data.frame(Cells = rownames(mat_H), mat_H)
  mat_H$Cells = factor(mat_H$Cells, levels = mat_H$Cells)
  
  
  mat_H <-
    mat_H %>% melt(
      id.vars = "Cells",
      variable.name = "GEPs",
      value.name = "Expression"
    )
  
  p <-
    ggplot(mat_H, aes(x = Cells, y = GEPs)) + geom_tile(aes(fill = Expression)) +
    scale_fill_gradient(
      name = "Other cells",
      low = "#F0FFFF",
      high = RColorBrewer::brewer.pal(9, "Blues")[8],
      limit = c(min(mat_H$Expression), 1),
      space = "Lab",
      guide = "colourbar"
    ) +
    new_scale("fill") +
    geom_tile(aes(fill = Expression), data = mat_H[which(mat_H$Cells %in% cell_ids),]) +
    scale_fill_gradient(
      name = "Cells of interest",
      low = "#F0FFFF",
      high = RColorBrewer::brewer.pal(9, "Reds")[8],
      limit = c(min(mat_H$Expression), 1),
      space = "Lab",
      guide = "colourbar"
    ) +
    xlab("Cells") +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title = element_text(size = 32, face = "bold", color = "black"),
      axis.ticks.length = unit(.20, "cm"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        size = 24,
        face = "bold",
        colour = "black"
      ),
      axis.text.x = element_blank(),
      legend.title = element_text(
        size = 24,
        face = "bold",
        colour = "black"
      ),
      legend.key.size = unit(3, "line"),
      legend.text = element_text(size = 24, face = "bold")
    )
  ggsave(
    filename = str_c(temp_dir, "Coefficient_heatmap_", figureName, ".png"),
    plot = p,
    width = 44,
    height = 28,
    dpi = 300
  )
  
  
}
