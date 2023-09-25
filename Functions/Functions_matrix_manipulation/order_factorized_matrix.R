# To scale and order the basis and coefficient matrix
order_factorized_matrix <- function(mat) {
  df <- mat
  df = as.data.frame(t(apply(df, 1, function(x) {
    x / sum(x)
  })))
  df_list = list()
  
  for (i in(c(1:ncol(df)))) {
    temp = df[apply(df, 1, function(x) {
      which.max(x)
    }) == i,]
    temp = temp[order(temp[, i], decreasing = TRUE), ]
    df_list[[i]] <- temp
  }
  df <- do.call(rbind.data.frame, df_list)
  return(df)
}
