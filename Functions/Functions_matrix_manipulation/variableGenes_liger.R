# This function takes normalized count :: For Liger, it's relative count of cells. Cell count is divided by the cell's total count
# This function takes a seurat object or a matrix (for now, I will be using only seurat object as input of the function)

variableGenes_liger <- function(dataObject, 
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
  
  # tc_per_cell (total read count per cell) :: it calculates the total raw read count of cells
  tc_per_cell <- colSums(dataObject@assays$RNA@counts)
  
  # Each gene's mean expression level (across all cells)
  gene_expr_mean <- rowMeansFast(dataObject@assays$RNA@data)
  
  # Each gene's expression variance (across all cells)
  gene_expr_var <- rowVarsFast(dataObject@assays$RNA@data, gene_expr_mean)
  
  # Assign gene names to the calculated mean and variance values
  names(gene_expr_mean) <- names(gene_expr_var) <- rownames(dataObject@assays$RNA@data)
  
  # Find an explanation for this nolan constant :::: Must be done
  nolan_constant <- mean((1 / tc_per_cell))
  
  # Why alpha.thresh was corrected :::: Must be done
  alphathresh.corrected <- alpha.thresh / nrow(dataObject)
  
  # Why gene mean upper bound was introduced? And how exactly is it calculated? :::: Must be done
  genemeanupper <- gene_expr_mean + qnorm(1 - alphathresh.corrected / 2) *
    sqrt(gene_expr_mean * nolan_constant / ncol(dataObject))
  
  # Why gene mean lower bound was introduced? And how exactly is it calculated? :::: Must be done
  basegenelower <- log10(gene_expr_mean * nolan_constant)
  
  # This function identifies genes that are above certain threshold. Filter genes to keep the desired number of genes specified by the user.
  num_varGenes <- function(x, num.genes.des){
    # This function returns the difference between the desired number of genes and
    # the number actually obtained when thresholded on x
    y <- length(which(gene_expr_var / nolan_constant > genemeanupper &
                        log10(gene_expr_var) > basegenelower + x))
    return(abs(num.genes.des - y))
  }
  
  # Optimize the value of "var.thresh" to keep the desired number of genes 
  optimized <- optimize(num_varGenes, c(0, 1.5), tol = tol, num.genes.des = num.genes)
  
  # Assign the new "var.thresh" value based on the optimization step performed above
  var.thresh <- optimized$minimum
  
  if (optimized$objective > 1) {
    warning(paste0("Returned number of genes for dataset ", 1, " differs from requested by ",
                   optimized$objective, ". Lower tol or alpha.thresh for better results."))}
  
  # Get the gene names that were selected as highly variable features by Liger.
  # Explain how is it done :::: Must be done
  genes.new <- names(gene_expr_var)[which(gene_expr_var / nolan_constant > genemeanupper &
                                            log10(gene_expr_var) > basegenelower + var.thresh)]
  
  
  return(genes.new)
}
