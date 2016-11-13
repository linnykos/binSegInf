#' Generate polyhedra matrix from bsFs object
#' 
#' Forms both Gamma matrix and u vector
#'
#' @param obj bsFs object
#' @param y numeric vector to represent data
#' @param ... not used
#'
#' @return An object of class polyhedra
#' @export
form_polyhedra.bsFs <- function(obj, y, ...){
  if(.get_startEnd(obj$tree$name)[2] != length(y)) stop("obj does not match y")
  isValid(obj)
  
  n <- length(y)
  numSteps <- obj$numSteps
  comp.lis <- .list_comparison(obj)
  gamma.row.lis <- vector("list", numSteps)
  
  hash_nodes <- hash::hash()
  
  for(i in 1:numSteps){
    losing.mat <- comp.lis[[i]]$losing
    
    gamma.row.lis[[i]] <- .gammaRows_from_comparisons(comp.lis[[i]]$winning,
      losing.mat, y)
  }
  
  mat <- do.call(rbind, gamma.row.lis)
  polyhedra(gamma = mat, u = rep(0, nrow(mat)))
}

.gammaRows_from_comparisons <- function(vec, mat, y){
  stopifnot(length(vec) == 3, ncol(mat) == 3)

  n <- length(y)
  win.contrast <- .cusum_contrast_full(vec[1], vec[2], vec[3], n)
  lose.contrast <- t(apply(mat, 1, function(x){
    .cusum_contrast_full(x[1], x[2], x[3], n)
  }))
  
  sign.win <- as.numeric(sign(win.contrast %*% y))
  
  # add inequalities to compare winning split to all other splits
  res <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win, 
    rep(1, nrow(lose.contrast)))
  res2 <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win, 
    -rep(1, nrow(lose.contrast)))
  
  # add inequalities to compare splits to 0 (ensure correct sign)
  rbind(res, res2, sign.win * win.contrast)
}

.vector_matrix_signedDiff <- function(vec, mat, sign.vec, sign.mat){
  stopifnot(!is.matrix(vec), is.numeric(vec), is.matrix(mat), is.numeric(mat))
  stopifnot(length(vec) == ncol(mat))
  stopifnot(length(sign.vec) == 1, length(sign.mat) == nrow(mat))
  stopifnot(all(c(sign.vec, sign.mat) %in% c(-1,0,1)))
  
  t(sign.vec * vec - t(sign.mat * mat))
}