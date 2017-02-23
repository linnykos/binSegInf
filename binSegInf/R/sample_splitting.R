#' Sample splitting 
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
  res <- jumps(method(y[idx], ...))
  
  #readjust the changepoints
  res <- (res-1)*2 + 1
  
  #average the points on the second half of the data
  structure(list(jumps = res, n = n, method = deparse(substitute(method))), class = "ss")
}

#' Get jumps from bsFs objects
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
  
  if(jump > n-3) warning("The detected jump is too close to the boundary")
  
  #figure out the indices to use in the contrast
  jumpSorted.vec <- c(1, jumps(obj, T), n)
  idx <- which(jumpSorted.vec == jump)[1]
  if(idx == 1) idx <- 2
  start <- jumpSorted.vec[idx-1]; split <- jump; end <- jumpSorted.vec[idx+1]
  
  pos_idx <- seq(jump+3, ceiling((end-1)/2)*2, by = 2)
  neg_idx <- seq(floor((start+1)/2)*2, jump-1, by = 2)
  
  v <- rep(0, n)
  v[pos_idx] <- 1/length(pos_idx); v[neg_idx] <- -1/length(neg_idx)
  
  v
}

#' P-values for post-selection inference
#'
#' @param y numeric vector
#' @param contrast contrast numeric vector
#'
#' @return a numeric p-value between 0 and 1
#' @export
pvalue_ss <- function(y, contrast){
  pos_idx <- which(contrast > 0); neg_idx <- which(contrast < 0)
  n_pos <- length(pos_idx); n_neg <- length(neg_idx)
  pos_mean <- mean(y[pos_idx]); neg_mean <- mean(y[neg_idx])
  pos_sd <- stats::sd(y[pos_idx]); neg_sd <- stats::sd(y[neg_idx])
  
  if(is.na(pos_sd)) pos_sd <- 0
  if(is.na(neg_sd)) neg_sd <- 0
  
  z <- (pos_mean - neg_mean)/(sqrt(pos_sd/n_pos + neg_sd/n_neg))
  2*stats::pnorm(-abs(z))
}