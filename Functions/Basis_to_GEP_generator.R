# Get GEPs from the generated Basis file - containing grouping of the genes
Basis_to_GEP_generator <- function(input_data, 
                          Factor_ID = "50", 
                          store_GEPs = TRUE,
                          return_GEPs = FALSE,
                          generate_ortho = FALSE,
                          ortho_file = NULL,
                          conversion_ID = NULL,
                          GEP_IDs = NULL,
                          store_dir = NULL, 
                          store_folder = "Factor_Ks_GEPs"){
  
  # Creating necessary storing space to store the results
  
  # Output directory
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Output folder
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # Create separate directory to store the ortho GEPs
  if (generate_ortho == TRUE) {
    
    # Output folder for specific factorization rank
    if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Factor_K_", Factor_ID, "_GEPs_Ortho"))){
      dir.create(str_c(store_dir, "/", store_folder, "/", "Factor_K_", Factor_ID, "_GEPs_Ortho"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Variable containing the path to the output folder
    temp_dir = str_c(store_dir, "/", store_folder, "/", "Factor_K_", Factor_ID, "_GEPs_Ortho", "/")
  }
  
  # Create separate directory to store the ortho GEPs
  if (generate_ortho == FALSE) {
    # Output folder for specific factorization rank
    if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "Factor_K_", Factor_ID, "_GEPs"))){
      dir.create(str_c(store_dir, "/", store_folder, "/", "Factor_K_", Factor_ID, "_GEPs"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Variable containing the path to the output folder
    temp_dir = str_c(store_dir, "/", store_folder, "/", "Factor_K_", Factor_ID, "_GEPs", "/")
  }
  
  
  # Column name in the basis file
  Factor_col = colnames(input_data)[grep(pattern = Factor_ID, x = colnames(input_data))]
  
  # Cluster levels of the genes / GEP ids for the factorization rank
  Factor_GEPs = levels(input_data[, Factor_col])
  
  # Creating list containing genes in each GEP
  GEPs_list = sapply(Factor_GEPs, function(x) {rownames(input_data)[which(input_data[, Factor_col, drop = FALSE] == x)]}, simplify = FALSE, USE.NAMES = TRUE)
  
  # Setting names of the items in the list
  names(GEPs_list) <- str_c("GEP_", names(GEPs_list))
  
  # Check if user called for particular GEP IDs
  if (!missing(GEP_IDs)){
    GEPs_list = GEPs_list[sapply(GEP_IDs, function(x) {grep(pattern = str_c("^", x, "$", sep = ""), parse_number(names(GEPs_list)))}, simplify = "vector")]
  }
  
  # Check if user specified to store the outputs
  if (store_GEPs == TRUE) {
    
    # Check if the user wants orthologues as output
    if (generate_ortho == TRUE) {
      
      # Check if the ortho file is provided
      if (missing(ortho_file)) {
        print("No orthologues table found. Please provide an orthologues table with gene IDs to look for in the gene_ID column.")
      }
      
      else {
        # Transform the gene IDs to orthologues IDs
        GEPs_list <- lapply(GEPs_list, function(gene_list) {
          return(ortho_file[ortho_file[, "gene_ID"] %in% gene_list, conversion_ID])
        })
        
        # Setting names of the items in the list
        names(GEPs_list) <- str_c("GEP_", names(GEPs_list))
        
        # Write the GEPs into text files
        invisible(lapply(names(GEPs_list), function (i) {fileGenerator(GEPs_list[[i]], fileName = str_c(temp_dir, i, "_ortho.txt"))}))
        
        # If user wants to return the output file alongside storing the GEPs
        if (return_GEPs == TRUE) {
          
          # return the GEP list
          return(GEPs_list)
        }
      }
    }
    
    else {
      # Write the GEPs into text files without returning any output in the console
      invisible(lapply(names(GEPs_list), function (i) {fileGenerator(GEPs_list[[i]], fileName = str_c(temp_dir, i, ".txt"))}))
      
      if (return_GEPs == TRUE) {
        # return the GEP list
        return(GEPs_list)
      }
      
    }
    
  }
  
  # Only returning the output without writing any files
  else {
    
    # Check if the user wants orthologues as output
    if (generate_ortho == TRUE) {
      
      # Check if the ortho file is provided
      if (missing(ortho_file)) {
        print("No orthologues table found. Please provide an orthologues table with gene IDs to look for in the gene_ID column.")
      }
      
      else {
        
        # Transform the gene IDs to orthologues IDs
        GEPs_list <- lapply(GEPs_list, function(gene_list) {
          return(ortho_file[ortho_file[, "gene_ID"] %in% gene_list, conversion_ID])
        })
        
        # return the GEP list
        return(GEPs_list)
        }
      }
    
    else {
      
      # return the GEP list
      return(GEPs_list)
    }
  }
}
