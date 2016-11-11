polyhedra <- function(gamma, u){
  structure(list(gamma = gamma, u = u), class = "polyhedra")
}

isValid.polyhedra <- function(obj){
  if(!is.numeric(obj$gamma) | !is.matrix(obj$gamma)) 
    stop("gamma is not a numeric matrix")
  if(!is.numeric(obj$u)) stop("u is not a numeric vector")
  if(nrow(obj$gamma) != length(obj$u)) stop("nrow(gamma) does not match length(u)")
  
  TRUE
}