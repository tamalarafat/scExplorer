# Transform the dataset to ortho table :: replace gene ids of a given species with ortho genes of target species
prepare_ortho_data <- function(input_data, ortho_data, ortho_column_name_of_gene_ids, ortho_column_name_to_assign){
  # get the gene ids in the data
  temp_gene_ids = rownames(input_data)
  
  # get the gene ids in the ortho table
  temp_ortho_ids = as.character(ortho_data[[ortho_column_name_of_gene_ids]])
  
  # This step is optional. All the hirsuta gene ids that are present in the orthologues table should be present in Cardamine genome. But anyway, if there is anything wrong with the naming of the genes, this step will take care of it.
  temp_ortho_ids <- intersect(temp_gene_ids, temp_ortho_ids)
  
  # Subsetting the data with the orthogenes only
  input_data <- input_data[temp_ortho_ids, ]
  
  # Create a vector of column names that will be dropped later from the dataframe
  drop_columns <- c("CHID", ortho_column_name_to_assign) ## Columns to drop
  
  # Converting the data into a data.frame so I can do whatever I want.
  temp_input_data <- as.data.frame(input_data)
  
  # Create a column with the rownames of the data. Later we will use this information to merge the data tables (scRNAseq - ortho table)
  temp_input_data$CHID <- rownames(temp_input_data)
  
  # Merge the two tables
  temp_input_data <- merge(x = temp_input_data, y = ortho_data, by.x = "CHID", by.y = ortho_column_name_of_gene_ids, all.x = T)
  
  # assign the row names with the desired ids
  rownames(temp_input_data) <- temp_input_data[ , ortho_column_name_to_assign]
  
  # convert the table to a sparse matrix
  temp_input_data <- as(as.matrix(temp_input_data[ ,!(names(temp_input_data) %in% drop_columns)]), "sparseMatrix")
  
  return(temp_input_data)
}
