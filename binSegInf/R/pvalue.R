pvalue <- function(y, polyhedra, contrast, sigma = 1, value = 0){
  terms <- .compute_truncGaus_terms(y, polyhedra, contrast, sigma)
  
  sapply(value, function(x){.truncated_gauss_cdf(terms$term, mu = x, 
    sigma = terms$sigma, a = terms$a, b = terms$b)})
}

.compute_truncGaus_terms <- function(y, polyhedra, contrast, sigma){
  z <- as.numeric(contrast %*% y)
  
  if(z < 0) contrast <- -contrast; z <- as.numeric(contrast %*% y)
  
  vv <- contrast %*% contrast
  sd <- as.numeric(sigma*sqrt(vv))
  
  rho <- as.numeric(polyhedra$gamma %*% contrast) / vv
  vec <- as.numeric((polyhedra$u - polyhedra$gamma %*% y + rho * z)/rho)
  
  if(any(rho > 0)) vlo <- max(vec[rho > 0]) else vlo <- -Inf
  if(any(rho < 0)) vup <- min(vec[rho < 0]) else vup <- Inf
  
  list(term = z, sigma = sd, a = vlo, b = vup)
}

.truncated_gauss_cdf <- function(value, mu, sigma, a, b, tol = 1e-4){
  sapply(value, function(x){
    if(x <= a) { 
      0
    } else if(x >= b){
      1
    } else {
      a <- Rmpfr::mpfr((a-mu)/sigma, precBits = 10)
      b <- Rmpfr::mpfr((b-mu)/sigma, precBits = 10)
      z <- Rmpfr::mpfr((value-mu)/sigma, precBits = 10)
      as.numeric((Rmpfr::pnorm(b) - Rmpfr::pnorm(z))/(Rmpfr::pnorm(b) - Rmpfr::pnorm(a)))
    }
  })
}