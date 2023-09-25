# This function performs GO annotation for any provided set of gene ids or a data.frame or a list of data.frames
markers_GO_annotation <- function(marker_set,
                                  species_gene_set,
                                  output_file_name = "GO_annotation_of_genes",
                                  store_dir = NULL, 
                                  store_folder = "RES_0_",
                                  GO_results_folder_name = "Cluster_specific_conserved_marker") {
  
  # marker_set; Input marker file, can be a character vector, or a dataframe, or a list of multiple dataframes
  # species_gene_set; All gene ids of the species of interest
  # output_file_name; if the input is a character vector or a dataframe, user can specify the name of the output file
  # store_dir; path of the directory to store the files
  # store_folder; the specified folder will be created in the specified directory with the provided name to store the files; for a particular resolution, all the results for the resolution parameter will be stored here
  # GO_results_folder_name; A folder will be created with the user specified name to store the GO annotation results for a particular run
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder, "/", GO_results_folder_name))){
    dir.create(str_c(store_dir, "/", store_folder, "/", GO_results_folder_name), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/", store_folder, "/", GO_results_folder_name, "/")
  
  # Clean the query cache of biomaRt
  biomaRt::biomartCacheClear()
  
  # Creating an empty list to store the gene ids of interest
  temp_list = list()
  
  # Let's check the input file if it's a character vector or a data.frame or a list of dataframes
  # Optional :: We can make one additional criteria for the list items, where the list items can be multiple character vector of genes 
  
  # If the input is a character vector, we convert it to a list
  if (class(marker_set) == "character") {
    temp_list[[1]] = marker_set
    names(temp_list)[1] = output_file_name
  } 
  # If the input is a data.frame, we get the row names where the gene ids are embedded and convert it to a list
  else if (class(marker_set) == "data.frame") {
    temp_marker_list = rownames(marker_set)
    temp_list[[1]] = temp_marker_list
    names(temp_list)[1] = output_file_name
  } 
  # If the input is a list containing multiple data.frames, we get the row names for each data.frame and stroe it in the list
  else if (class(marker_set) == "list") {
    if (is.null(names(marker_set))){
      names(marker_set) <- str_c("module_", c(1:length(marker_set)))
    }
    
    temp_names = names(marker_set)
    
    for (i in c(1:length(temp_names))){
      
      if (is.character(marker_set[[temp_names[i]]])) {
        temp_marker_list = marker_set[[temp_names[i]]]
        temp_list[[i]] = temp_marker_list
        names(temp_list)[i] = temp_names[i]
      }
      else {
        temp_marker_list = rownames(marker_set[[temp_names[i]]])
        temp_list[[i]] = temp_marker_list
        names(temp_list)[i] = temp_names[i]
      }
    }
    
  } else {
    print("The input file can be a character vector containing gene IDs, or a dataframe with row names containing gene IDs, or a list containing multiple dataframes with row names containing gene IDs.")
  }
  
  
  # information retrieval ::  TAIR / Ensembl / NCBI Entrez gene ids for the marker ids
  attributes_to_retrieve = c("tair_symbol", "uniprotswissprot", "entrezgene_id")
  
  # For all the genes in the species of interest, query the Ensembl API, and retrieve ids from biomartr 
  species_genes_annotated <- biomartr::biomart(genes = species_gene_set,
                                               mart       = "plants_mart",                 
                                               dataset    = "athaliana_eg_gene",           
                                               attributes = attributes_to_retrieve,        
                                               filters    =  "ensembl_gene_id")
  
  # Get the item names of the list to generate folder with the name
  temp_item_names = names(temp_list)
  
  for (j in c(1:length(temp_item_names))){
    
    # Let's create a directory to store the GO annotation files
    dir.create(str_c(temp_dir, temp_item_names[j]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    
    # Assigning the directory path to a local variable
    temp_go_dir = str_c(temp_dir, temp_item_names[j], "/")
    
    # Let's get the gene ids from the list
    temp_marker_genes = temp_list[[temp_item_names[j]]]
    
    ## Biological processes
    
    # Lets retreive the attributes for the marker set
    GeneSet_annotated <- biomartr::biomart(genes = temp_marker_genes,
                                           mart       = "plants_mart",                 
                                           dataset    = "athaliana_eg_gene",           
                                           attributes = attributes_to_retrieve,        
                                           filters =  "ensembl_gene_id" )
    
    # performing the over representation analysis (ORA) for Gene Ontology class - Biological processes
    GeneSet_ORA_BP <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                               universe = species_genes_annotated$entrezgene_id,
                               OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                               keyType = "ENTREZID", 
                               ont = "BP",
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05,
                               readable = TRUE,
                               pool = FALSE)
    
    tryCatch({
      # Lets save the GO annotation results for BP class
      GeneSet_all_BP = GeneSet_ORA_BP@result
      
      # All GO terms (significant, non-significant)
      write.csv(GeneSet_all_BP, file = str_c(temp_go_dir, temp_item_names[j], "_markers_All_BP.csv"))
      
      # The Gene Ontology classification is very redundant; meaning that parental terms overlap a lot with their related child terms.
      # The clusterProfiler package comes with a dedicated function called "simplify" to solve this issue.
      # We can also manually extract the significant terms based on setting a p.adjust threshold of < 0.05
      
      GeneSet_ORA_BP_simplified <- clusterProfiler::simplify(GeneSet_ORA_BP)
      
      GeneSet_qurated_BP = GeneSet_ORA_BP_simplified@result
      
      # Simplified annotation files contain only the significant GO terms; q value < 0.05
      write.csv(GeneSet_qurated_BP, file = str_c(temp_go_dir, temp_item_names[j], "_markers_Significant_BP.csv"))
    }, 
    error = function(e){str_c("No biological processe annotation was found for ", temp_item_names[j], " markers.")}
    )
    
    # Molecular functions
    
    # performing the ORA for Gene Ontology class - Molecular function
    GeneSet_ORA_MF <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                               universe = species_genes_annotated$entrezgene_id,
                               OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                               keyType = "ENTREZID", 
                               ont = "MF",
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05,
                               readable = TRUE,
                               pool = FALSE)
    
    tryCatch({
      # Lets save the GO annotation results for MF class
      GeneSet_all_MF = GeneSet_ORA_MF@result
      
      write.csv(GeneSet_all_MF, file = str_c(temp_go_dir, temp_item_names[j], "_markers_All_MF.csv"))
      
      GeneSet_ORA_MF_simplified <- clusterProfiler::simplify(GeneSet_ORA_MF)
      
      GeneSet_qurated_MF = GeneSet_ORA_MF_simplified@result
      
      write.csv(GeneSet_qurated_MF, file = str_c(temp_go_dir, temp_item_names[j], "_markers_Significant_MF.csv"))
    }, 
    error = function(e){str_c("No molecular function annotation was found for ", temp_item_names[j], " markers.")}
    )
    
    # Cellular components
    
    # performing the ORA for Gene Ontology class - Cellular component
    GeneSet_ORA_CC <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                               universe = species_genes_annotated$entrezgene_id,
                               OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                               keyType = "ENTREZID", 
                               ont = "CC",
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05,
                               readable = TRUE,
                               pool = FALSE)
    
    tryCatch({
      # Lets save the GO annotation results for CC class
      GeneSet_all_CC = GeneSet_ORA_CC@result
      
      write.csv(GeneSet_all_CC, file = str_c(temp_go_dir, temp_item_names[j], "_markers_All_CC.csv"))
      
      GeneSet_ORA_CC_simplified <- clusterProfiler::simplify(GeneSet_ORA_CC)
      
      GeneSet_qurated_CC = GeneSet_ORA_CC_simplified@result
      
      write.csv(GeneSet_qurated_CC, file = str_c(temp_go_dir, temp_item_names[j], "_markers_Significant_CC.csv"))
    }, 
    error = function(e){str_c("No cellular component was found for ", temp_item_names[j], " markers.")}
    )
    
    ### Visualization of the enricchment analysis results
    
    # Figure 1 :: Generate dotplot for all GO terms - BP
    if (dim(GeneSet_ORA_BP)[1] < 1){
      print("No module found")
    } 
    else {
      p1 <- dotplot(GeneSet_ORA_BP, showCategory = 20) + 
        theme(axis.title.x = element_text(size = 16, face = "bold"), 
              axis.ticks.length = unit(.30, "cm"), 
              axis.text = element_text(size = 72, face = "bold"),
              legend.key = element_rect(size = 28),
              legend.text = element_text(size = 16, face = "bold"),
              legend.title = element_text(size = 16, face = "bold"))
      
      ggsave(filename = str_c(temp_go_dir, "GO_annotation_ALL_BP_",temp_item_names[j], "_markers.png"), plot = p1, width = 12, height = 14, dpi = 300)
    }
    
    
    # Figure 2 :: Generate dotplot for significant GO terms - BP
    if (dim(GeneSet_ORA_BP_simplified)[1] < 1){
      print("No module found")
    } 
    else {
      p2 <- dotplot(GeneSet_ORA_BP_simplified, showCategory = 20) + 
        theme(axis.title.x = element_text(size = 16, face = "bold"), 
              axis.ticks.length = unit(.30, "cm"), 
              axis.text = element_text(size = 72, face = "bold"),
              legend.key = element_rect(size = 28),
              legend.text = element_text(size = 16, face = "bold"),
              legend.title = element_text(size = 16, face = "bold"))
      
      ggsave(filename = str_c(temp_go_dir, "GO_annotation_Significant_BP_",temp_item_names[j], "_markers.png"), plot = p2, width = 12, height = 14, dpi = 300)
    }
    
    # Figure 3 :: Generate network of GO terms using all GO terms - BP
    if (dim(GeneSet_ORA_BP)[1] > 0){
      # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. This makes it easier to identify functional modules.
      GeneSet_ORA_BP <- pairwise_termsim(GeneSet_ORA_BP, method = "JC")
      
      if (dim(GeneSet_ORA_BP@termsim)[1] == 1){
        print("No module found")
      } 
      else {
        tryCatch({
          p3 <- emapplot(GeneSet_ORA_BP, color = "qvalue")  + theme_bw() + 
            theme(axis.title = element_blank(), 
                  axis.ticks.length = unit(.30, "cm"), 
                  axis.text = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, face = "bold"))
          
          ggsave(filename = str_c(temp_go_dir, "GO_modules_ALL_BP_",temp_item_names[j], "_markers.png"), plot = p3, width = 12, height = 14, dpi = 300)
        }, 
        error = function(e){"No module found"}
        )
      }
    }
    
    else {
      print("No module found")
    }
    
    # Figure 4 :: Generate network of GO terms using significant GO terms - BP
    if (dim(GeneSet_ORA_BP_simplified)[1] > 0){
      # For the qurated or simplifeid BP
      # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. This makes it easier to identify functional modules.
      GeneSet_ORA_BP_simplified <- pairwise_termsim(GeneSet_ORA_BP_simplified, method = "JC")
      
      if (dim(GeneSet_ORA_BP_simplified@termsim)[1] == 1){
        print("No module found")
      } 
      else {
        tryCatch({
          p4 <- emapplot(GeneSet_ORA_BP_simplified, color = "qvalue") + theme_bw() + 
            theme(axis.title = element_blank(), 
                  axis.ticks.length = unit(.30, "cm"), 
                  axis.text = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, face = "bold"))
          
          ggsave(filename = str_c(temp_go_dir, "GO_modules_Significant_BP_",temp_item_names[j], "_markers.png"), plot = p4, width = 12, height = 14, dpi = 300)
        }, 
        error = function(e){"No module found"}
        )
      }
    }
    
    else {
      print("No module found")
    }
    
  }
}