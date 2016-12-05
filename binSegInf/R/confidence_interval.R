#' Confidence interval from post-selection inference
#'
#' @param y numeric vector
#' @param polyhedra polyhedra object
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param alpha numeric between 0 and 1 with default of 0.05
#' @param gridsize numeric to denote how fine of a grid to invert the hypothesis
#' test
#' @param alternative string of either "one.sided" or "two.sided" for the 
#' alternative. If one.sided, the alternative means the test statistic is positive.
#' @param precBits the number of bits used to compute the p value
#'
#' @return a vector of two numbers, the lower and upper end of the confidence interval
#' @export
confidence_interval <- function(y, polyhedra, contrast, sigma = 1, alpha = 0.05,
  gridsize = 250, alternative = c("two.sided", "one.sided"), precBits = 10){
  alternative <- match.arg(alternative, c("two.sided", "one.sided"))
  
  diff <- max(y) - min(y)
  seq.val <- seq(-2*diff, 2*diff, length.out = gridsize)
  pvalue <- pvalue(y, polyhedra, contrast, sigma, null_mean = seq.val,
    precBits = precBits)
  
  if(alternative == "two.sided"){
    idx <- c(.select_index(pvalue, alpha/2, T), .select_index(pvalue, 1-alpha/2, F))
    c(seq.val[idx[1]], seq.val[idx[2]])
  } else {
    idx <- .select_index(pvalue, alpha, T)
    c(seq.val[idx], Inf)
  }
}
  
.select_index <- function(vec, alpha, lower = T){
  idx <- ifelse(lower, min(which(vec >= alpha)), max(which(vec <= alpha)))
  if(length(idx) == 0 | is.na(idx) | is.infinite(idx)){
    warning("Numeric precision suspected to be too low")
    if(lower) return(1) else return(length(vec))
  }

  if(lower & vec[idx] > alpha & idx > 1) idx <- idx - 1
  if(!lower & vec[idx] < alpha & idx < length(vec)) idx <- idx + 1
  
  idx
}