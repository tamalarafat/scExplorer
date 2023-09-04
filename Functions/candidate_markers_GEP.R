# Get GEPs from the generated Basis file - containing grouping of the genes
candidate_markers_GEP <- function(seurat_object, 
                                  GEP_genes = NULL, 
                                  GEP_IDs = NULL, 
                                  cluster_ID,
                                  DE_test = "wilcox", 
                                  reduction_name = "inmf",
                                  find_candidates = 10,
                                  GEP_source_name = NULL,
                                  combine_GEPs = FALSE,
                                  store_outputs = TRUE,
                                  store_dir = NULL, 
                                  store_folder = "GEP_candidates_of_clusters"){
  
  # Creating necessary storing space to store the results
  
  if (store_outputs == TRUE) {
    # Output directory
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    # Output folder
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Output folder for the target cluster
    if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "cluster_", cluster_ID))){
      dir.create(str_c(store_dir, "/", store_folder, "/", "cluster_", cluster_ID), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Assigning the output directory path to a variable
    temp_dir = str_c(store_dir, "/", store_folder, "/", "cluster_", cluster_ID, "/")
  }
  
  
  
  # Get the genes of the GEP of interest. User can specify the genes as
  # 1 - Without specifying any GEP IDs or set of genes. If not specified, the characteristic GEP of the cluster will be picked to get the genes
  # 2 - File path including the text file. Only .txt is allowed
  # 3 - A vector containing the gene IDs
  # 4 - GEP ID or IDs
  
  # This part will take the argument input and get the genes of the GEP or GEPs and assign them into a list that will be used later to iterate over the GEP IDs.
  
  # 1- If GEP IDs or GEP genes is not provided - get the characteristic GEP of the target cluster
  if (missing(GEP_IDs) & missing(GEP_genes)) {
    
    # If GEP ID or IDs are not provided, use the characteristic GEP of the cell cluster to identify candidate markers
    defining_GEP = characteristic_GEP_of_cells(seurat_object = seurat_object, target_ID = cluster_ID, store_output = FALSE)
    
    # Get the ID of the characteristic GEP
    GEP_ID = as.character(defining_GEP[which.max(defining_GEP$cell_count), "GEP"])
    
    # Get the genes belonging to the identified GEP
    GEP_list = get_GEP_genes(seurat_object = seurat_object, reduction_name = reduction_name, GEP_IDs = GEP_ID)
    
  } 
  
  # 2 & 3- If file path or a vector of genes provided
  else if (!missing(GEP_genes)) {
    
    # Check if the input GEP_genes is a path with the file name. Only .txt file is allowed
    if (file_test("-f", GEP_genes)) {
      
      # Read the txt file and add colnames
      GEP_genes = read.delim(GEP_genes, header = FALSE, col.names = "gene_ID")
      
      # Get the genes and create a list with the gene IDs
      GEP_list = list(GEP_genes = GEP_genes$gene_ID)
    } 
    
    # Check if the input GEP_genes is a vector containing the gene IDs
    else {
      
      # create a list with the gene IDs
      GEP_list = list(GEP_genes = GEP_genes)
    }
  } 
  
  # 4- If GEP ID or IDs are provided 
  else if (!missing(GEP_IDs)){
    
    # Get the genes of the provided GEP ids. This function returns a list
    GEP_list = get_GEP_genes(seurat_object = seurat_object, reduction_name = reduction_name, GEP_IDs = GEP_IDs)
  }
  
  # extract the basis matrix from the Seurat object
  basis = as.data.frame(seurat_object@reductions[[reduction_name]]@feature.loadings)
  
  # Add column names to the basis matrix
  colnames(basis) = str_c("GEP_", parse_number(colnames(basis)))
  
  # Order the matrix with all geps
  basis = order_factorized_matrix(basis)
  
  # Create an empty list to store the candidate markers for each GEP ID
  candidates_list = list()
  
  # Iterate over the list items. List items hold the genes in the GEPs
  for (i in c(1:length(GEP_list))) {
    
    # Get the genes in each item / GEP
    GEPgenes = GEP_list[[i]]
    
    # Differential expression test to perform. By default, uses the "wilcox" test provided by the Seurat tool
    
    # Run differential expression test for the genes in the desired GEP.
    cluster_GEP_markers = FindMarkers(seurat_object, ident.1 = cluster_ID, features = GEPgenes, test.use = DE_test, only.pos = TRUE)
    
    # Add the gene names to the DE file
    cluster_GEP_markers$gene_ID = rownames(cluster_GEP_markers)
    
    # Subset the basis matrix with the GEP genes
    df = basis[GEPgenes, ]
    
      
    # If missing GEP IDs, extract the GEP ID from the basis matrix
    if(missing(GEP_IDs)){
      GEP_ID = names(which.max(table(apply(df, 1, which.max))))
    } 
    
    # Else take the provided GEP ID
    else {
      GEP_ID = GEP_IDs[i]
    }
  
    # Subset the basis matrix with the constituting genes of the GEP
    df = df[, GEP_ID, drop = FALSE]
    
    # Set the column names to coefficient
    colnames(df)[1] <- "coefficient"
    
    # Order the rows or genes based on the coefficient value
    df = df[order(df$coefficient, decreasing = TRUE), , drop = FALSE]
    
    # Add rank information - ranks are based on the coefficient. Largest value top rank, lowest value bottom rank
    df$rank = c(1:nrow(df))
    
    # Set the rownames with the gene IDs
    df$gene_ID = rownames(df)
    
    
    # Create a data.frame by merging the coefficient and differential expression analysis results - candidates data.frame
    candidates_df = merge(df, cluster_GEP_markers, by = "gene_ID")
    
    # Order the data.frame according to the pct.2 criteria where lowest pct corresponds to the top row / order
    candidates_df = candidates_df[order(candidates_df$pct.2, decreasing = FALSE), ]
    
    # Subset the data.frame with the pct.2 criteria. Keep only those rows with less than 0.1 (10%) pct.2 detection
    candidates_df = candidates_df[candidates_df$pct.2 < 0.1, ]
    
    # if the selection criteria does not work, we pick the top 10 genes
    
    # Check if the data.frame is empty. No genes have less than 0.1 pct.2 detection
    if (nrow(candidates_df) == 0) {
      paste("No genes")
      } 
    
    # Check if number of row is more than one and less than or equal to user requested number of candidate genes. If Yes, pick all the rows
    # This part ensures that all the genes with less than 10% (0.1) pct.2 detection is picked as candidates
    else if (nrow(candidates_df) <= find_candidates & nrow(candidates_df) != 0){
        candidates_df = candidates_df
      } 
    
    # If we have more than requested number of candidate genes with <0.1 pct.2 detection, we consider the average fold change criteria
    else {
        if (nrow(candidates_df[candidates_df$avg_log2FC >= 1,]) < find_candidates & nrow(candidates_df[candidates_df$avg_log2FC >= 1,]) > 0){
          subset1 = candidates_df[candidates_df$avg_log2FC >= 1,]
          subset2 = candidates_df[candidates_df$avg_log2FC < 1,]
          subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
          
          if (nrow(subset2) <= (find_candidates - nrow(subset1))) {
            candidates_df = rbind.data.frame(subset1, subset2)
          } 
          
          else {
            subset2 = head(subset2, (find_candidates - nrow(subset1)))
            candidates_df = rbind.data.frame(subset1, subset2)
          }
        }
      
        else if (nrow(candidates_df[candidates_df$avg_log2FC >= 1,]) == 0){
          candidates_df = candidates_df[order(candidates_df$avg_log2FC, decreasing = TRUE), ]
          candidates_df = candidates_df[c(1:find_candidates), ]
        }
      
        else {
          candidates_df = candidates_df[candidates_df$avg_log2FC >= 1, ]
          candidates_df = candidates_df[c(1:find_candidates), ]
        }
      }
    
    
    # If the condition above finds less than desired number of candidate genes, that this part will be executed
    if (nrow(candidates_df) < find_candidates){
      
      candidate_subset = merge(df, cluster_GEP_markers, by = "gene_ID")
      
      candidate_subset = candidate_subset[order(candidate_subset$pct.2, decreasing = FALSE), ]
      
      candidate_subset = candidate_subset[!candidate_subset$gene_ID %in% rownames(candidates_df), ]
      
      if (nrow(candidate_subset) <= (find_candidates - nrow(candidates_df))) {
        candidates_df = rbind.data.frame(candidates_df, candidate_subset)
      } 
      
      else {
        candidate_subset = head(candidate_subset, (find_candidates - nrow(candidates_df)))
        candidates_df = rbind.data.frame(candidates_df, candidate_subset)
      }
      
    }
    
    rownames(candidates_df) <- candidates_df$gene_ID
    
    if (!missing(GEP_source_name)){
      
      candidates_df$source = GEP_source_name
    }
    
    else {
      candidates_df$source = str_c("GEP_", GEP_ID)
    }
    
    candidates_list[[i]] <- candidates_df
    
    names(candidates_list)[i] <- GEP_ID
  }
  
  if (combine_GEPs == FALSE) {
    
    if (store_outputs == TRUE) {
      
      for (i in length(candidates_list)) {
        write.csv(candidates_list[[i]], file = str_c(temp_dir, "Cluster_", cluster_ID, "_", names(candidates_list)[i],"_candidates.csv"), row.names = FALSE)
      } 
    }
    else {
      return(candidates_list)
    }
  }
  
  else {
    
    candidates_list = do.call(rbind.data.frame, candidates_list)
    
    if (store_outputs == TRUE) {
      
      write.csv(candidates_list, file = str_c(temp_dir, "Cluster_", cluster_ID, "_", "GEP_", paste(GEP_IDs, collapse = "_"),"_candidates.csv"), row.names = FALSE)
    }
    
    else {
      return(candidates_list)
    }
  }
}
