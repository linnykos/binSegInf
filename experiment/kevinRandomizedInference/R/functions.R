#' Estimate the largest directional difference in means
#'
#' @param y data vector
#'
#' @return an index (positive integer)
#' @export
estimate_kevin <- function(y){
NULL
}

#' Compute polyhedraon of \code{estimate_kevin}
#'
#' @param y data vector
#' @param i index where \code{estimate_kevin} estimated
#'
#' @return a matrix to be the Gamma matrix
#' @export
polyhedron_kevin <- function(y, i){
  NULL
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
