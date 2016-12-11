#' P-values for post-selection inference
#'
#' @param y numeric vector
#' @param polyhedra polyhedra object
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param null_mean the null-hypothesis mean to test against
#' @param alternative string of either "one.sided" or "two.sided" for the 
#' alternative. If one.sided, the alternative means the test statistic is positive.
#' @param precBits the number of bits used to compute the p value
#'
#' @return a numeric p-value between 0 and 1
#' @export
pvalue <- function(y, polyhedra, contrast, sigma = 1, null_mean = 0,
  alternative = c("one.sided", "two.sided"), precBits = 10){
  alternative <- match.arg(alternative, c("one.sided", "two.sided"))
  
  terms <- .compute_truncGaus_terms(y, polyhedra, contrast, sigma)
  
  res <- sapply(null_mean, function(x){.truncated_gauss_cdf(terms$term, mu = x, 
    sigma = terms$sigma, a = terms$a, b = terms$b, precBits = precBits)})
  
  if(alternative == "one.sided") res else 2*min(res, 1-res)
}

.compute_truncGaus_terms <- function(y, polyhedra, contrast, sigma){
  z <- as.numeric(contrast %*% y)
  
  vv <- contrast %*% contrast
  sd <- as.numeric(sigma*sqrt(vv))
  
  rho <- as.numeric(polyhedra$gamma %*% contrast) / vv
  vec <- as.numeric((polyhedra$u - polyhedra$gamma %*% y + rho * z)/rho)
  
  if(any(rho > 0)) vlo <- max(vec[rho > 0]) else vlo <- -Inf
  if(any(rho < 0)) vup <- min(vec[rho < 0]) else vup <- Inf
  
  list(term = z, sigma = sd, a = vlo, b = vup)
}

.truncated_gauss_cdf <- function(value, mu, sigma, a, b, tol_zero = 1e-5, 
  tol_prec = 1e-2, precBits = 10){
  if(b < a) stop("b must be greater or equal to a")
  
  sapply(value, function(x){
    if(x <= a) return(0)
    if(x >= b) return(1)
   
    a_scaled <- (a-mu)/sigma; b_scaled <- (b-mu)/sigma
    z_scaled <- (value-mu)/sigma
    denom <- stats::pnorm(b_scaled) - stats::pnorm(a_scaled)
    numerator <- stats::pnorm(b_scaled) - stats::pnorm(z_scaled)
    
    if(denom > tol_prec & numerator > tol_prec & numerator/denom > tol_prec 
      & numerator/denom < 1-tol_prec){
      return(numerator/denom)
    } else {
      a_scaled <- Rmpfr::mpfr((a-mu)/sigma, precBits = precBits)
      b_scaled <- Rmpfr::mpfr((b-mu)/sigma, precBits = precBits)
      z_scaled <- Rmpfr::mpfr((value-mu)/sigma, precBits = precBits)  
  
      denom <- Rmpfr::pnorm(b_scaled) - Rmpfr::pnorm(a_scaled)
      if(denom < tol_zero) denom <- tol_zero
      as.numeric((Rmpfr::pnorm(b_scaled) - Rmpfr::pnorm(z_scaled))/denom)
    }
  })
}
