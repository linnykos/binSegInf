polyhedra <- function(gamma, u){
  structure(list(gamma = gamma, u = u), class = "polyhedra")
}

isValid.polyhedra <- function(obj){
  if(nrow(obj$gamma) != length(obj$u)) stop("nrow(gamma) does not match length(u)")
  
  TRUE
}