# Convert the coefficient matrix to a s4 object
transform_coef_to_S4 <- function(nmfObject){
  # The coefficient matrix
  coef_mat <- t(nmfObject@fit@H)
  
  setClass("DimReduc", representation(cell.embeddings = "matrix", feature.loadings = "matrix", assay.used = "character", global = "logical", key = "character"))
  nnmf <- new("DimReduc", cell.embeddings = coef_mat, assay.used = 'RNA')
  nnmf
}
