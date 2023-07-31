# This function takes normalized count :: For Liger, it's relative count of cells. Cell count is divided by the cell's total count
# This function takes a seurat object or a matrix (for now, I will be using only seurat object as input of the function)

generate_mean_var_liger <- function(dataObject, 
                                    var.thresh = 0.1,
                                    num.genes = 3000,
                                    alpha.thresh = 0.99,
                                    tol =  0.0001
){
  
  # var.thresh :: Variance threshold. Main threshold used to identify variable genes. Genes with expression variance greater than threshold (relative to mean) are selected.
  # (higher threshold -> fewer selected genes). Accepts single value or vector with separate var.thresh for each dataset. (default 0.1)
  
  # num.genes :: Number of genes to find for the dataset. Optimizes the value of var.thresh for each dataset to get this number of genes.
  
  # alpha.thresh :: Controls upper bound for expected mean gene expression (lower threshold -> higher upper bound). (default 0.99)
  
  # tol :: Tolerance to use for optimization if num.genes values passed in (default 0.0001).
  
  # Calculate mean of the genes
  rowMeansFast <- function(x) {
    .Call('_rliger_rowMeansFast', PACKAGE = 'rliger', x)
  }
  
  # Calculate variance in respect to the mean expression of the genes
  rowVarsFast <- function(x, means) {
    .Call('_rliger_rowVarsFast', PACKAGE = 'rliger', x, means)
  }
  
  # Each gene's mean expression level (across all cells)
  gene_expr_mean <- rowMeansFast(dataObject@assays$RNA@data)
  
  # Each gene's expression variance (across all cells)
  gene_expr_var <- rowVarsFast(dataObject@assays$RNA@data, gene_expr_mean)
  
  # Assign gene names to the calculated mean and variance values
  names(gene_expr_mean) <- names(gene_expr_var) <- rownames(dataObject@assays$RNA@data)
  
  df_mean_var <- data.frame(avgExpression = log10(gene_expr_mean), geneVariance = log10(gene_expr_var))
  
  return(df_mean_var)
  
}
