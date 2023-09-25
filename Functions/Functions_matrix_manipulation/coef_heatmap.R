### Coefficient heatmap generator
coef_heatmap <- function(matList, figureName) {
  
  if (class(matList) != "list") {
    mat_H = data.frame(Cells = rownames(matList), matList)
    mat_H$Cells = factor(mat_H$Cells, levels = mat_H$Cells)
    
    mat_H <-
      mat_H %>% melt(
        id.vars = "Cells",
        variable.name = "GEPs",
        value.name = "Expression"
      )
    
    p <-
      ggplot(mat_H, aes(x = Cells, y = GEPs, fill = Expression)) + geom_tile() +
      scale_fill_gradient(
        name = "Expression level",
        low = RColorBrewer::brewer.pal(9, "Blues")[1],
        high = RColorBrewer::brewer.pal(9, "Blues")[8],
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
        axis.text.y = element_text(size = 24, face = "bold", colour = "black"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 24, face = "bold", colour = "black"),
        legend.key.size = unit(3, "line"),
        legend.text = element_text(size = 24, face = "bold")
      )
    
    ggsave(filename = figureName, plot = p, width = 44, height = 28, dpi = 300)
  }
  
  else {
    Datasets <- names(matList)
    
    for (d in c(1:length(Datasets))) {
      mat_H <- matList[[Datasets[d]]]
      mat_H = data.frame(Cells = rownames(mat_H), mat_H)
      mat_H$Cells = factor(mat_H$Cells, levels = mat_H$Cells)
      
      mat_H <-
        mat_H %>% melt(
          id.vars = "Cells",
          variable.name = "GEPs",
          value.name = "Expression"
        )
      
      p <-
        ggplot(mat_H, aes(x = Cells, y = GEPs, fill = Expression)) + geom_tile() +
        scale_fill_gradient(
          name = "Expression level",
          # changes legend title
          low = RColorBrewer::brewer.pal(9, "Blues")[1],
          high = RColorBrewer::brewer.pal(9, "Blues")[8],
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
          # Background of the entire plot
          axis.title = element_text(
            size = 26,
            face = "bold",
            color = "black"
          ),
          axis.ticks.length = unit(.20, "cm"),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
          axis.text.x = element_blank(),
          legend.title = element_text(size = 20, face = "bold", colour = "black"),
          legend.key.size = unit(2, "line"),
          legend.text = element_text(size = 20, face = "bold")
        )
      
      ggsave(filename = str_c("Coefficient_heatmap_", Datasets[d], "_", figureName, ".png"), plot = p, width = 44, height = 28, dpi = 300)
    }
  }
}
