### Generates heatmap of the dataset-specific basis matrix
basis_heatmap <- function(matList, figureName) {
  
  if (class(matList) != "list") {
    mat_W = data.frame(Genes = rownames(matList), matList)
    mat_W$Genes <- factor(mat_W$Genes, levels = mat_W$Genes)
    mat_W <-
      mat_W %>% melt(
        id.vars = "Genes",
        variable.name = "GEPs",
        value.name = "Expression"
      )
    
    p <- ggplot(mat_W, aes(x = GEPs, y = Genes, fill = Expression)) + geom_tile() +
      scale_fill_gradient(
        name = "Expression level",
        low = RColorBrewer::brewer.pal(9, "Blues")[1],
        high = RColorBrewer::brewer.pal(9, "Blues")[8],
        limit = c(min(mat_W$Expression), 1),
        space = "Lab",
        guide = "colourbar"
      ) +
      xlab("GEPs") +
      theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 32, face = "bold", color = "black"),
        axis.ticks.length = unit(.20, "cm"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 24, face = "bold", colour = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 24, face = "bold", colour = "black"),
        legend.key.size = unit(3, "line"),
        legend.text = element_text(size = 24, face = "bold")
      )
    
    ggsave(filename = figureName, plot = p, width = 32, height = 32, dpi = 300)
  }
  
  else {
    Datasets <- names(matList)
    
    for (d in c(1:length(Datasets))) {
      mat_W <- matList[[Datasets[d]]]
      mat_W = data.frame(Genes = rownames(mat_W), mat_W)
      mat_W$Genes <- factor(mat_W$Genes, levels = mat_W$Genes)
      mat_W <-
        mat_W %>% melt(
          id.vars = "Genes",
          variable.name = "GEPs",
          value.name = "Expression"
        )
      
      p <-
        ggplot(mat_W, aes(x = GEPs, y = Genes, fill = Expression)) + geom_tile() +
        scale_fill_gradient(
          name = "Expression level",
          # changes legend title
          low = RColorBrewer::brewer.pal(9, "Blues")[1],
          high = RColorBrewer::brewer.pal(9, "Blues")[8],
          limit = c(min(mat_W$Expression), 1),
          space = "Lab",
          guide = "colourbar"
        ) +
        xlab("GEPs") +
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
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(
            size = 20,
            face = "bold",
            colour = "black",
            angle = 90,
            vjust = 0.5
          ),
          axis.text.y = element_blank(),
          legend.title = element_text(
            size = 20,
            face = "bold",
            colour = "black"
          ),
          legend.key.size = unit(2, "line"),
          legend.text = element_text(size = 20, face = "bold")
        )
      ggsave(filename = str_c("Basis_heatmap_", Datasets[d], "_", figureName, ".png"), plot = p, width = 32, height = 32, dpi = 300)
    }
  }
}