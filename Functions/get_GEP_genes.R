get_GEP_genes <- function (seurat_object,
                           reduction_name = "inmf",
                           GEP_IDs = NULL) {
  
  # Extract the basis matrix
  temp_mat <- as.data.frame(seurat_object@reductions[[reduction_name]]@feature.loadings)
  
  # Add columns to store the grouping of the genes to different GEPs.
  temp_mat = temp_mat %>% 
    mutate(belonging_GEP = colnames(temp_mat)[apply(temp_mat, 1, which.max)], GEPs = parse_number(colnames(temp_mat)[apply(temp_mat, 1, which.max)])) %>% 
    dplyr::select(c(belonging_GEP, GEPs)) %>% mutate(GEPs = as.factor(GEPs))
  
  # Create an empty list to store the genes of each GEP
  GEPs_list = list()
  
  GEP_labels = str_sort(as.character(levels(temp_mat$GEPs)), numeric = TRUE)
  
  # iterate over the GEP ids
  for (i in c(1:length(GEP_labels))) {
    
    # get the gene IDs from the data.frame
    gene_names = rownames(temp_mat[temp_mat$GEPs == i, ])
    
    # Add the gene names in the list
    GEPs_list[[i]] <- gene_names
    
    # Add names to the list item
    names(GEPs_list)[i] <- str_c("GEP_", i)
    
  }

  # Check if user called for particular GEP IDs
  if (!missing(GEP_IDs)){
    GEPs_list = GEPs_list[sapply(GEP_IDs, function(x) {grep(pattern = str_c("^", x, "$", sep = ""), parse_number(names(GEPs_list)))}, simplify = "vector")]
    }
  
  # return the GEP list
  return (GEPs_list)
  
}
