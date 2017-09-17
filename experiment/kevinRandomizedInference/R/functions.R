#' Estimate the largest directional difference in means
#'
#' @param y data vector
#'
#' @return an index (positive integer)
#' @export
estimate_kevin <- function(y){
  n <- length(y)
  vec <- sapply(1:(n-1), function(i){
    mean(y[1:i]) - mean(y[(i+1):n])
  })
  which.max(vec)
}

#' Compute polyhedraon of \code{estimate_kevin}
#'
#' @param y data vector
#' @param i index where \code{estimate_kevin} estimated
#'
#' @return a matrix to be the Gamma matrix
#' @export
polyhedron_kevin <- function(y, i){
  n <- length(y)

  winning_contrast <- .polyhedron_vector_generator(i, n)

  idx <- c(1:(n-1))[-i]
  mat <- sapply(idx, function(x){
    losing_contrast <- .polyhedron_vector_generator(x, n)
    winning_contrast - losing_contrast
  })

  t(mat)
}

#' Compute pvalue from polyhedron
#'
#' @param mat Gamma matrix
#' @param y data vector
#' @param sigma2 variance of the generative process of y
#' @param contrast vector
#'
#' @return a list of the pvalue, the numerator of the pvalue, and the denominator of the pvalue
#' @export
poly.pval_kevin <- function(mat, y, sigma2, contrast){
  NULL
}

#' Sampler for randomized inference
#'
#' @param y data vector
#' @param sigma02 variance of the randomized noise being added
#' @param sigma2 variance of the generative process of y
#' @param num_trials number of trials for the sampler
#' @param contrast vector
#'
#' @return pvalue
#' @export
sampler_kevin <- function(y, sigma02, sigma2, num_trials, contrast){
  NULL
}

##########

.polyhedron_vector_generator <- function(i, n){
  vec <- rep(0, n)
  vec[1:i] <- 1/i
  vec[(i+1):n] <- -1/(n-i)
  vec
}
