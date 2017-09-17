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
  terms <- .compute_truncGaus_terms(y, mat, contrast, sigma2)

  res <- .truncated_gauss_cdf(terms$term, sigma = terms$sigma, a = terms$a,
                              b = terms$b)
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


.compute_truncGaus_terms <- function(y, mat, contrast, sigma2){
  z <- as.numeric(contrast %*% y)

  vv <- contrast %*% contrast
  sd <- as.numeric(sigma2*sqrt(vv))

  rho <- as.numeric(mat %*% contrast) / vv
  vec <- as.numeric((rep(0, nrow(mat)) - mat %*% y + rho * z)/rho)

  if(any(rho > 0)) vlo <- max(vec[rho > 0]) else vlo <- -Inf
  if(any(rho < 0)) vup <- min(vec[rho < 0]) else vup <- Inf

  list(term = z, sigma = sd, a = vlo, b = vup)
}

#' Title
#'
#' The mean that is being tested is always assumed to be 0.
#'
#' @param value
#' @param sigma
#' @param a
#' @param b
#' @param tol_prec
#'
#' @return
#' @export
#'
#' @examples
.truncated_gauss_cdf <- function(value, sigma2, a, b, tol_prec = 1e-2){
  if(b < a) stop("b must be greater or equal to a")

  if(val <= a) return(0)
  if(val >= b) return(1)

  a_scaled <- a/sigma2; b_scaled <- b/sigma2
  z_scaled <- value/sigma2
  denom <- stats::pnorm(b_scaled) - stats::pnorm(a_scaled)
  numerator <- stats::pnorm(b_scaled) - stats::pnorm(z_scaled)

  val <- numerator/denom

  #fix any NAs first
  issue <- any(is.na(val) | is.nan(val))
  if(any(issue)) val <- .truncated_gauss_cdf_Rmpfr(value, sigma2, a, b, 10)

  #fix any source of possible imprecision
  issue <- any(denom < tol_prec) | any(numerator < tol_prec) |
    any(val < tol_prec) | any(val > 1-tol_prec)

  if(any(issue) & !is.na(precBits)) val <- .truncated_gauss_cdf_Rmpfr(value, sigma2, a, b, 10)

  val
}

.truncated_gauss_cdf_Rmpfr <- function(value, sigma2, a, b, tol_zero = 1e-5,
                                       precBits = 10){

  a_scaled <- Rmpfr::mpfr(a/sigma2, precBits = precBits)
  b_scaled <- Rmpfr::mpfr(b/sigma2, precBits = precBits)
  z_scaled <- Rmpfr::mpfr(value/sigma2, precBits = precBits)

  denom <- Rmpfr::pnorm(b_scaled) - Rmpfr::pnorm(a_scaled)
  if(denom < tol_zero) denom <- tol_zero
  as.numeric(Rmpfr::mpfr((Rmpfr::pnorm(b_scaled) - Rmpfr::pnorm(z_scaled))/denom, precBits = precBits))
}

