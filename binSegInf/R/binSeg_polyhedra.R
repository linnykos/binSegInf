form_polyhedra.bsFs <- function(obj, y, ...){
  isValid(obj)
  
  n <- .extract_startEnd(obj$tree$name)[2]
  numSteps <- obj$numSteps
  comp.lis <- .list_comparison(obj)
  gamma.row.lis <- vector("list", numSteps)
  
  for(i in 1:numSteps){
    gamma.row.lis[[i]] <- .gammaRows_from_comparisons(comp.lis[[i]]$winning,
      comp.lis[[i]]$losing, n, y)
  }
  
  mat <- do.call(rbind, gamma.row.lis)
  polyhedra(gamma = mat, u = rep(0, nrow(mat)))
}

.gammaRows_from_comparisons <- function(vec, mat, n, y){
  stopifnot(length(vec) == 3, ncol(mat) == 3)

  winning.contrast <- .cusum_contrast_full(vec[1], vec[2], vec[3], n)
  losing.contrast <- t(apply(mat, 1, function(x){
    .cusum_contrast_full(x[1], x[2], x[3], n)
  }))
  
  sign.winning <- as.numeric(sign(winning.contrast %*% y))
  signs.losing <- as.numeric(sign(losing.contrast %*% y))
  
  # add inequalities to compare winning split to all other splits
  res <- .vector_matrix_signedDiff(winning.contrast, losing.contrast, sign.winning,
    signs.losing)
  
  # add inequalities to compare splits to 0 (ensure correct sign)
  rbind(res, sign.winning * winning.contrast, signs.losing * losing.contrast)
}

.vector_matrix_signedDiff <- function(vec, mat, sign.vec, sign.mat){
  stopifnot(!is.matrix(vec), is.numeric(vec), is.matrix(mat), is.numeric(mat))
  stopifnot(length(vec) == ncol(mat))
  stopifnot(length(sign.vec) == 1, length(sign.mat) == nrow(mat))
  stopifnot(all(c(sign.vec, sign.mat) %in% c(-1,0,1)))
  
  t(sign.vec * vec - t(sign.mat * mat))
}