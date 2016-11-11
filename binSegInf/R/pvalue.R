pvalue <- function(y, polyhedra, contrast, value = 0, sigma = 1){
  z <- as.numeric(contrast %*% y)
  
  if(z > 0) contrast <- -contrast; z <- as.numeric(contrast %*% y)
  
  vv <- contrast %*% contrast
  sd <- as.numeric(sigma*sqrt(vv))
  
  rho <- as.numeric(polyhedra$gamma %*% contrast) / vv
  vec <- as.numeric((polyhedra$u - polyhedra$gamma %*% y + rho * z)/rho)
  
  if(any(rho > 0)) vlo <- max(vec[rho > 0]) else vlo <- -Inf
  if(any(rho < 0)) vup <- min(vec[rho < 0]) else vup <- Inf
  
  sapply(value, function(x){.truncated_gauss_cdf(z, mean = x, 
    sigma = sd, a = vlo, b = vup)})
}

.truncated_gauss_cdf <- function(value, mean, sigma, a, b, tol = 1e-4){
  sapply(value, function(x){
    if(x <= a) { 
      0
    } else if(x >= b){
      1
    } else {
      denom <- stats::pnorm(b, mean, sigma) - stats::pnorm(a, mean, sigma)
      if(denom < tol) denom <- tol
      (stats::pnorm(x, mean, sigma) - stats::pnorm(a, mean, sigma))/denom
    }
  })
}