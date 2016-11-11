#' Generate polyhedra
#'
#' @param gamma numeric matrix
#' @param u numeric vector
#'
#' @return object of class polyhedra
#' @export
polyhedra <- function(gamma, u){
  structure(list(gamma = gamma, u = u), class = "polyhedra")
}

#' Check if polyhedra is valid
#'
#' @param obj object of class polyhedra
#'
#' @return TRUE if valid
#' @export
isValid.polyhedra <- function(obj){
  if(!is.numeric(obj$gamma) | !is.matrix(obj$gamma)) 
    stop("gamma is not a numeric matrix")
  if(!is.numeric(obj$u)) stop("u is not a numeric vector")
  if(nrow(obj$gamma) != length(obj$u)) stop("nrow(gamma) does not match length(u)")
  
  TRUE
}