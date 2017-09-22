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

  mat <- t(mat)
  list(gamma = mat, u = rep(0, nrow(mat)))
}

#' Compute pvalue from polyhedron
#'
#' @param y data vector
#' @param poly polyhedron object
#' @param contrast vector
#' @param sigma sd of the generative process of y
#'
#' @return a list of the pvalue, the numerator of the pvalue, and the denominator of the pvalue
#' @export
poly.pval_kevin <- function(y, poly, contrast, sigma){
  terms <- .compute_truncGaus_terms(y, poly, contrast, sigma)

  res <- .truncated_gauss_cdf(terms$term, sigma = terms$sigma, a = terms$a,
                              b = terms$b)
}

#' Sampler for randomized inference
#'
#' @param y data vector
#' @param sigma0 sd of the randomized noise being added
#' @param sigma sd of the generative process of y
#' @param num_trials number of trials for the sampler
#' @param contrast vector
#'
#' @return pvalue
#' @export
sampler_kevin <- function(y, noise, sigma0, sigma, num_trials, contrast){
  stopifnot(length(y) == length(noise))
  unique_seed <- sum(abs(y))

  n <- length(y)
  i <- estimate_kevin(y + noise)
  poly <- polyhedron_kevin(y + noise, i)

  vec <- sapply(1:num_trials, function(x){
    set.seed(x*10*unique_seed)
    new_noise <- stats::rnorm(n, sd = sigma0)
    poly$u <- poly$gamma %*% new_noise
    res <- poly.pval_kevin(y, poly, contrast, sigma)

    c(res$pvalue, res$denominator)
  })

  as.numeric(vec[1,]%*%vec[2,]/sum(vec[2,]))
}

##########

.polyhedron_vector_generator <- function(i, n){
  vec <- rep(0, n)
  vec[1:i] <- 1/i
  vec[(i+1):n] <- -1/(n-i)
  vec
}

.compute_truncGaus_terms <- function(y, poly, contrast, sigma){
  z <- as.numeric(contrast %*% y)

  vv <- contrast %*% contrast
  sd <- as.numeric(sigma*sqrt(vv))

  rho <- as.numeric(poly$gamma %*% contrast) / vv
  vec <- as.numeric(poly$u - poly$gamma %*% y + rho * z)/rho

  if(any(rho > 0)) vlo <- max(vec[rho > 0]) else vlo <- -Inf
  if(any(rho < 0)) vup <- min(vec[rho < 0]) else vup <- Inf

  list(term = z, sigma = sd, a = vlo, b = vup)
}

#' Computing pvalue with vlo and vhigh
#'
#' The mean that is being tested is always assumed to be 0.
#'
#' @param value the observed value to compute the density at
#' @param sigma sd of the truncated Gaussian
#' @param a left (lower) truncation of Gaussian
#' @param b right (upper) truncation of Gaussian
#' @param tol_prec require precision
#'
#' @return a pvalue, its numerator and its denominator
#' @export
.truncated_gauss_cdf <- function(value, sigma, a, b, tol_prec = 1e-5){
  if(b < a) stop("b must be greater or equal to a")

  a_scaled <- a/sigma; b_scaled <- b/sigma
  z_scaled <- value/sigma
  denom <- stats::pnorm(b_scaled) - stats::pnorm(a_scaled)

  if(value <= a){
    numerator <- 0
  } else if (value >= b){
    numerator <- 1
  } else {
    numerator <- stats::pnorm(b_scaled) - stats::pnorm(z_scaled)
  }

  val <- numerator/denom

  #fix any source of possible imprecision
  issue <- is.na(val) | is.nan(val) | denom < tol_prec | numerator < tol_prec |
    val < tol_prec | val > 1-tol_prec

  if(issue) {
    res <- .truncated_gauss_cdf_Rmpfr(value, sigma, a, b, 10)
    val <- res$val; numerator <- res$numerator; denom <- res$denom
  }

  list(pvalue = val, numerator = numerator, denominator = denom)
}

.truncated_gauss_cdf_Rmpfr <- function(value, sigma, a, b, tol_zero = 1e-8,
                                       precBits = 50){

  a_scaled <- Rmpfr::mpfr(a/sigma, precBits = precBits)
  b_scaled <- Rmpfr::mpfr(b/sigma, precBits = precBits)
  z_scaled <- Rmpfr::mpfr(value/sigma, precBits = precBits)

  denom <- Rmpfr::pnorm(b_scaled) - Rmpfr::pnorm(a_scaled)
  if(denom < tol_zero) denom <- tol_zero
  numerator <- as.numeric(Rmpfr::pnorm(b_scaled) - Rmpfr::pnorm(z_scaled))
  val <- as.numeric(Rmpfr::mpfr(numerator/denom, precBits = precBits))

  list(val = val, numerator = numerator, denom = denom)
}

