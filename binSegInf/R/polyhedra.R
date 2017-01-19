#' Generate polyhedra
#'
#' @param obj numeric matrix to be the gamma matrix
#' @param u numeric vector
#' @param ... not used
#'
#' @return object of class polyhedra
#' @export
polyhedra.matrix <- function(obj, u, ...){
  structure(list(gamma = obj, u = u), class = "polyhedra")
}

#' Check if polyhedra is valid
#'
#' @param obj object of class polyhedra
#'
#' @return TRUE if valid
#' @export
is_valid.polyhedra <- function(obj){
  if(!is.numeric(obj$gamma) | !is.matrix(obj$gamma)) 
    stop("gamma is not a numeric matrix")
  if(!is.numeric(obj$u)) stop("u is not a numeric vector")
  if(nrow(obj$gamma) != length(obj$u)) stop("nrow(gamma) does not match length(u)")
  
  TRUE
}