# Author - Yasir 
###
# Loading all the packages
###

libs = c(
  "Seurat",
  "clustree",
  "rliger",
  "cowplot",
  "ggplot2",
  "scExplorer",
  "ggpubr",
  "ggthemes",
  "ggsci",
  "scales",
  "openxlsx",
  "tidyr",
  "Matrix",
  "dplyr",
  "stringr",
  "readr",
  "wesanderson",
  "clusterProfiler",
  "slingshot",
  "viridis",
  "ggbeeswarm",
  "SingleCellExperiment",
  "stringi",
  "forcats",
  "biomartr",
  "enrichplot",
  "org.At.tair.db",
  "biomaRt",
  "NMF",
  "reshape2", 
  "ggnewscale",
  "CountClust",
  "ggrepel",
  "ggsankey"
)

invisible(lapply(libs, library, character.only = TRUE))

options(ggrepl.max.overlaps = Inf)

# rm(libs)
