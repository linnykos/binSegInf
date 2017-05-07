#pg 560 from https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxh008
circular_bin <- function(vec){
  n <- length(vec)
  cum_vec <- cumsum(vec)
  index_mat <- .form_index(n)
  
  vec <- .compute_cirBS(cum_vec, index_mat)
  max(abs(vec))
}

.compute_cirBS <- function(cum_vec, index_mat){
  stopifnot(is.matrix(index_mat), ncol(index_mat) == 2, min(index_mat) >= 1,
            max(index_mat) <= length(cum_vec))
  n <- length(cum_vec)
  
  sqrt(1/(index_mat[,2] - index_mat[,1]) + 1/(n - index_mat[,2] + index_mat[,1]))*
    ((cum_vec[index_mat[,2]] - cum_vec[index_mat[,1]])/(index_mat[,2] - index_mat[,1]) -
       (cum_vec[n] - cum_vec[index_mat[,2]] + cum_vec[index_mat[,1]])/(n - index_mat[,2] + index_mat[,1]))
}

.form_index <- function(n){
  index_1 <- rep(1:(n-1), times = (n-1):1)
  index_2 <- unlist(sapply(1:(n-1), function(x){(x+1):n}))
  cbind(index_1, index_2)
}