# Perform GO annotaation of the basis genes
Basis_genes_GO <- function(yourDF, ATgenes){
  
  # If biomart refuses to query Ensembl again, run this command:
  biomaRt::biomartCacheClear()
  
  # get TAIR/Ensembl/NCBI Entrez gene ids for the gene list
  attributes_to_retrieve = c("tair_symbol", "uniprotswissprot", "entrezgene_id")
  
  # For all the genes in the species of interest, query the Ensembl API, and retrieve ids from biomartr 
  all_arabidopsis_genes_annotated <- biomartr::biomart(genes = ATgenes,
                                                       mart       = "plants_mart",                 
                                                       dataset    = "athaliana_eg_gene",           
                                                       attributes = attributes_to_retrieve,        
                                                       filters    =  "ensembl_gene_id")
  
  # for compatibility with enrichGO universe genes in the universe need to be characters and not integers (Entrez gene id)
  all_arabidopsis_genes_annotated$entrezgene_id = as.character(all_arabidopsis_genes_annotated$entrezgene_id)
  
  # Go through each column of your dataframe
  # Each column contains information about on which genes belongs to which GEP / which set of coexpressed genes
  
  factorization_ID <- str_sub(colnames(yourDF), 3, -1) # Get the IDs/number of the factorization components, e.g. K = 3, 4, 5, 6, so on
  LF_ID <- colnames(yourDF) # Get the column names of yourDF / Basis components name
  
  for (i in c(1:ncol(yourDF))){
    
    # Let's create a directory to store the GO annotation results
    if (!dir.exists("GO_Annotation_of_Basis_Genes")){
      dir.create("GO_Annotation_of_Basis_Genes", showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Create a directory for the current iteration or factorization rank, e.g. for first i = 1, it will select the first column of the basis, K/LF = 3 (K = 3 components factorization)
    if (!dir.exists(str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i]))){
      dir.create(str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Let's store the factorization name in a temporary variable
    temp_LF_ID = LF_ID[i]
    
    # Go through each item or GEP of the selected column (i) of the basis - e.g. for K = 3, 3 GEPs, pick GEP 1, then 2, then 3
    for (j in c(1:length(levels(yourDF[, temp_LF_ID])))){
      
      # get the gene ids in the current GEP j (set of coexpressed genes)
      GeneSet <- rownames(yourDF[yourDF[, temp_LF_ID] == j, ])
      
      # Create a directory inside the previously created directory for the current GEP, j
      if (!dir.exists(str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j))){
        dir.create(str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j), showWarnings = TRUE, recursive = FALSE, mode = "0777")
      }
      
      # Save the set of genes in .txt file
      fileGenerator(GeneSet, fileName = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_genes.txt"))
      
      # Lets annotate the set of genes in the current GEP, j
      GeneSet_annotated <- biomartr::biomart(genes = GeneSet,
                                             mart       = "plants_mart",                 
                                             dataset    = "athaliana_eg_gene",           
                                             attributes = attributes_to_retrieve,        
                                             filters =  "ensembl_gene_id" )
      
      # performing the over representation analysis (ORA) for Gene Ontology class - Biological processes
      GeneSet_ORA_BP <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                                 universe = all_arabidopsis_genes_annotated$entrezgene_id,
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
        
        # All GO terms (significant, non-siignificant)
        write.csv(GeneSet_all_BP, file = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_BP_all.csv"))
        
        # The Gene Ontology classification is very redundant; meaning that parental terms overlap a lot with their related child terms.
        # The clusterProfiler package comes with a dedicated function called "simplify" to solve this issue.
        # We can also manually extract the significant terms based on setting a p.adjust threshold of < 0.05
        
        GeneSet_ORA_BP_simplified <- clusterProfiler::simplify(GeneSet_ORA_BP)
        
        GeneSet_qurated_BP = GeneSet_ORA_BP_simplified@result
        
        # Simplified annotation files contain only the significant GO terms; q value < 0.05
        write.csv(GeneSet_qurated_BP, file = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_significant_BP.csv"))
      }, 
      error = function(e){str_c("No biological processes was found for", " LF ", factorization_ID[i], " GEP ", j)}
      )
      
      # performing the ORA for Gene Ontology class - Molecular function
      GeneSet_ORA_MF <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                                 universe = all_arabidopsis_genes_annotated$entrezgene_id,
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
        
        write.csv(GeneSet_all_MF, file = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_MF_all.csv"))
        
        GeneSet_ORA_MF_simplified <- clusterProfiler::simplify(GeneSet_ORA_MF)
        
        GeneSet_qurated_MF = GeneSet_ORA_MF_simplified@result
        
        write.csv(GeneSet_qurated_MF, file = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_significant_MF.csv"))
      }, 
      error = function(e){str_c("No molecular functions was found for", " LF ", factorization_ID[i], " GEP ", j)}
      )
      
      # performing the ORA for Gene Ontology class - Cellular component
      GeneSet_ORA_CC <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                                 universe = all_arabidopsis_genes_annotated$entrezgene_id,
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
        
        write.csv(GeneSet_all_CC, file = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_CC_all.csv"))
        
        GeneSet_ORA_CC_simplified <- clusterProfiler::simplify(GeneSet_ORA_CC)
        
        GeneSet_qurated_CC = GeneSet_ORA_CC_simplified@result
        
        write.csv(GeneSet_qurated_CC, file = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_significant_CC.csv"))
      }, 
      error = function(e){str_c("No cellular component was found for", " LF ", factorization_ID[i], " GEP ", j)}
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
        
        ggsave(filename = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_all_BP_GO.png"), plot = p1, width = 12, height = 14, dpi = 300)
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
        
        ggsave(filename = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_sigificant_BP_GO.png"), plot = p2, width = 12, height = 14, dpi = 300)
      }
      
      # Figure 3 :: Generate network of GO terms using all GO terms - BP
      if (dim(GeneSet_ORA_BP)[1] > 0){
        # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. This makes it easier to identify functional modules.
        GeneSet_ORA_BP <- pairwise_termsim(GeneSet_ORA_BP, method = "JC")
        
        if (dim(GeneSet_ORA_BP@termsim)[1] == 1){
          print("No module found")
        } 
        else {
          p3 <- emapplot(GeneSet_ORA_BP, color = "qvalue")  + theme_bw() + 
            theme(axis.title = element_blank(), 
                  axis.ticks.length = unit(.30, "cm"), 
                  axis.text = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, face = "bold"))
          
          ggsave(filename = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_all_BP_GO_modules_network.png"), plot = p3, width = 12, height = 14, dpi = 300)
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
          p4 <- emapplot(GeneSet_ORA_BP_simplified, color = "qvalue") + theme_bw() + 
            theme(axis.title = element_blank(), 
                  axis.ticks.length = unit(.30, "cm"), 
                  axis.text = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, face = "bold"))
          
          ggsave(filename = str_c("GO_Annotation_of_Basis_Genes/Factorization_K_", factorization_ID[i], "/GEP_", j, "/GEP_", j, "_significant_BP_GO_modules_network.png"), plot = p4, width = 12, height = 14, dpi = 300)
        }
      }
      
      else {
        print("No module found")
      }
      
    }
  } 
}
