#' Sample splitting 
#'
#' Even points are for finding the changepoints. Odd points are used for testing.
#'
#' @param y data vector
#' @param method estimator on the first split
#' @param ... additional arguments
#'
#' @return ss object
#' @export
sample_splitting <- function(y, method, ...){
  n <- length(y)
  idx <- seq(2, n, by = 2)
  
  #estimate on half the data
  res <- jump_sign(method(y[idx], ...))
  
  #readjust the changepoints
  res$Jump <- res$Jump*2
  
  #average the points on the second half of the data
  structure(list(jumps = res$Jump, sign = res$Sign, n = n, 
                 method = deparse(substitute(method))), class = "ss")
}

#' Get jumps from ss objects
#' 
#' Enumerates the jumps. Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj bsFs object
#' @param sorted boolean
#' @param ... not used
#'
#' @return vector of jumps
#' @export
jumps.ss <- function(obj, sorted = F, ...){
  idx <- obj$jumps
  if(sorted) sort(idx) else idx
}

#' Generate contrast bector
#' 
#' Only does segment tests
#'
#' @param obj object
#' @param jump.idx index among the list of jumps to jump at
#' @param sorted boolean on whether or not jumps should be sorted
#'
#' @return numeric vector
#' @export
contrast_vector_ss <- function(obj, jump.idx, sorted = F){
  n <- obj$n
  jump.vec <- jumps(obj, sorted)
  jump <- jump.vec[jump.idx]
  jump_sign <- obj$sign[which(obj$jump == jump)]

  #figure out the indices to use in the contrast
  jumpSorted.vec <- c(1, jumps(obj, T), n)
  idx <- which(jumpSorted.vec == jump)[1]
  if(idx == 1) idx <- 2
  start <- jumpSorted.vec[idx-1]; split <- jump; end <- jumpSorted.vec[idx+1]
  
  if(jump + 3 < floor(end/2)*2-1){ 
    pos_idx <- seq(jump + 3, floor(end/2)*2-1, by = 2)
  } else {
    pos_idx <- min(jump + 3, floor((n+1)/2)*2-1)
  }
  
  if(floor(start/2)*2+3 < jump-1){
    neg_idx <- seq(floor(start/2)*2+3, jump-1, by = 2)
  } else {
    neg_idx <- max(jump - 1, 1)
  }
  
  v <- rep(0, n)
  v[pos_idx] <- 1/length(pos_idx); v[neg_idx] <- -1/length(neg_idx)
  
  attr(v, "sign") <- jump_sign
  v
}

#' P-values for sample splitting
#' 
#' This is limited to one-sided p-values currently, and the null-mean
#' is set to be 0.
#'
#' @param y numeric vector
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#'
#' @return a numeric p-value between 0 and 1
#' @export
pvalue_ss <- function(y, contrast, sigma = 1){
  pos_idx <- which(contrast > 0); neg_idx <- which(contrast < 0)
  n_pos <- length(pos_idx); n_neg <- length(neg_idx)
  pos_mean <- mean(y[pos_idx]); neg_mean <- mean(y[neg_idx])
  
  z <- (pos_mean - neg_mean)/(sqrt(sigma^2/n_pos + sigma^2/n_neg))
  
  if(attr(contrast, "sign") == -1) z <- -z
  
  1 - stats::pnorm(z)
}

#' Confidence interval for sample splitting
#' 
#' This is limited to two-sided confidence intervals currently
#'
#' @param y numeric vector
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param alpha numeric between 0 and 1 with default of 0.95. This is the
#' significance level, guaranteeing that alpha percentage
#' of the intervals will cover the true parameter.
#'
#' @return 
#' @export
confidence_interval_ss <- function(y, contrast, sigma = 1, alpha = 0.95){
  pos_idx <- which(contrast > 0); neg_idx <- which(contrast < 0)
  n_pos <- length(pos_idx); n_neg <- length(neg_idx)
  pos_mean <- mean(y[pos_idx]); neg_mean <- mean(y[neg_idx])
  
  interval_length <- abs(stats::qnorm((1-alpha)/2))*sqrt(sigma^2/n_pos + sigma^2/n_neg)
  c((pos_mean - neg_mean) - interval_length, (pos_mean - neg_mean) + interval_length)
}