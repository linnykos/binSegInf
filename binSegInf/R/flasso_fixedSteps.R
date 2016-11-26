flasso_fixedSteps <- function(y, numSteps, tol = 1e-7){
  if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
  if(any(duplicated(y))) stop("y must contain all unique values")

  #initialization
  n <- length(y)
  model.mat <- as.data.frame(matrix(NA, numSteps + 1, 3))
  colnames(model.mat) <- c("Index", "Sign", "Lambda")
  D <- .form_Dmatrix(n)

  for(steps in 1:(numSteps + 1)){
    idx <- .select_nonactive(n, model.mat$Index)
    a.vec <- .compute_fused_numerator(D, idx, y)
    b.vec <- .compute_fused_denominator(D, idx, model.mat[1:(steps-1),,drop = F])
    
    pos.ratio <- a.vec/(1+b.vec)
    pos.ratio[abs(1 + b.vec) < tol] <- 0
    neg.ratio <- a.vec/(-1+b.vec)
    neg.ratio[abs(-1 + b.vec) < tol] <- 0 

    if(max(pos.ratio) > max(neg.ratio)){
      model.mat[steps,] <- c(idx[which.max(pos.ratio)], 1, max(pos.ratio))
    } else {
      model.mat[steps,] <- c(idx[which.max(neg.ratio)], -1, max(neg.ratio))
    }
  }
  
  y.fit <- .refit_flasso(y, model.mat)

  structure(list(model = model.mat[1:numSteps, ], y.fit = y.fit, 
    numSteps = numSteps), class = "flFs")
}

jumps.flFs <- function(obj, sorted = F, ...){
  idx <- obj$model$Index
  if(sorted) sort(idx) else idx
}

jump_lambda.flFs <- function(obj, ...){
  obj$model$Lambda
}

summary.flFs <- function(object, ...){
  object$model
}

.refit_flasso <- function(y, model.mat){
  n <- length(y)
  D <- .form_Dmatrix(n)
  u.vec <- rep(NA, n-1)
  lambda <- min(model.mat$Lambda)
  u.vec[model.mat$Index] <- lambda*model.mat$Sign
  
  idx <- .select_nonactive(n, model.mat$Index)
  a.vec <- .compute_fused_numerator(D, idx, y)
  b.vec <- .compute_fused_denominator(D, idx, model.mat)
  
  nonactive.idx <- c(1:(n-1))[-model.mat$Index]
  u.vec[nonactive.idx] <- a.vec - lambda*b.vec
  
  as.numeric(y - t(.form_Dmatrix(n))%*%u.vec)
}

.form_Dmatrix <- function(n){
  t(sapply(1:(n-1), function(x){
    vec <- rep(0, n)
    vec[c(x, x+1)] <- c(-1,1)
    vec
  }))
}

.select_nonactive <- function(n, vec){
  val <- vec[!is.na(vec)]
  if(length(val) == 0) 1:(n-1) else c(1:(n-1))[-val]
}

.compute_fused_numerator <- function(D, idx, y){
  stopifnot(is.numeric(D), is.matrix(D), is.numeric(y))
  stopifnot(all(idx %% 1 == 0), !any(duplicated(idx)))
  stopifnot(min(idx) >= 1, max(idx) <= length(y) - 1)
  stopifnot(ncol(D) == length(y), nrow(D) == length(y) - 1)
  
  DDT <- D[idx,,drop = F]%*%t(D[idx,,drop = F])
  Dy <- D[idx,,drop = F]%*%y
  
  .svd_solve(DDT, Dy)
}

.compute_fused_denominator <- function(D, idx, model.mat){
  stopifnot(is.numeric(D), is.matrix(D))
  stopifnot(all(idx %% 1 == 0), !any(duplicated(idx)))
  stopifnot(min(idx) >= 1, max(idx) <= ncol(D) - 1)
  stopifnot(ncol(D) == nrow(D) + 1)
  
  if(length(idx) == nrow(D) || length(model.mat) == 0 || any(is.na(model.mat$Index))) 
    return(rep(0, nrow(D)))
  
  active.idx <- model.mat$Index; sign.vec <- model.mat$Sign
  DDT <- D[idx,,drop = F]%*%t(D[idx,,drop = F])
  DDTs <- D[idx,,drop = F] %*% t(D[active.idx,,drop = F]) %*% sign.vec
  
  .svd_solve(DDT, DDTs)
}

#solves Ax = b for A as a PSD matrix. Equivalently, (A.inv)b
.svd_solve <- function(A, b, tol = 1e-7){
  stopifnot(is.matrix(A), is.numeric(A), is.numeric(b))
  stopifnot(nrow(A) == length(b))
  
  s <- svd(A)
  d <- s$d
  stopifnot(all(d > -tol)) #ensure A is PSD
  bool <- (d > tol)
  d[bool] <- 1/d[bool]; d[!bool] <- 0
  
  as.numeric(s$v %*% (d * t(s$u) %*% b))
}