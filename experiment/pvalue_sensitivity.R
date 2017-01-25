pval_standard <- function(a, b, z){
  num <- stats::pnorm(as.numeric(b)) - stats::pnorm(as.numeric(z))
  denom <- stats::pnorm(as.numeric(b)) - stats::pnorm(as.numeric(a))
  num/denom
}

pval_Rmpfr <- function(a, b, z, precBits = 10){
  num <- Rmpfr::pnorm(as.numeric(b)) - Rmpfr::pnorm(as.numeric(z))
  denom <- Rmpfr::pnorm(as.numeric(b)) - Rmpfr::pnorm(as.numeric(a))
  as.numeric(Rmpfr::mpfr(num/denom, precBits = 10))
}

dif_vec <- seq(1e-3, 1, length.out = 100)
a_vec <- seq(0, 5, length.out = 100)

mat <- expand.grid(a_vec, dif_vec)
res <- apply(mat, 1, function(x){
  c(pval_standard(x[1], x[1]+10*x[2], x[1]+x[2]), pval_Rmpfr(x[1], x[1]+10*x[2], x[1]+x[2]))
})
res <- t(res)
dif_vec <- res[,1] - res[,2]
quantile(abs(dif_vec))
