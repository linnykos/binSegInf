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
  if(!is.numeric(obj$gamma) | !is.matrix(obj$gamma)) {
    print(head(obj$gamma))
    stop("gamma is not a numeric matrix")
  }
  if(!is.numeric(obj$u)) stop("u is not a numeric vector")
  if(nrow(obj$gamma) != length(obj$u)) stop("nrow(gamma) does not match length(u)")

  TRUE
}


#' Generic for functions that combine same-type things.
#' @param obj object
combine <- function(obj, ...) {UseMethod("combine")}

##' Combines several polyhedra to a single polyhedron
##' @param ... polyhedra objects to add
combine.polyhedra <- function(...){

    polyobjs = list(...)
    polyobjs = polyobjs[which(!sapply(polyobjs, is.null))]

    ## Combine separately and return
    newgamma = do.call(rbind, lapply(polyobjs, function(mypolyobj) mypolyobj$gamma))
    newu = as.numeric(do.call(c, lapply(polyobjs, function(mypolyobj) mypolyobj$u)))
    newpoly = structure(list(gamma = newgamma, u= newu), class = "polyhedra")

    stopifnot(is_valid.polyhedra(newpoly))
    return(newpoly)
}


print <- function(obj, ...) {UseMethod("print")}
## Prints polyhedra
print.polyhedra <- function(mypoly){
    if(all(is.na(mypoly$gamma[1,])) & nrow(mypoly$gamma)==1) print("Empty polyhedra object!")
    first.n = 10
    print(paste("Gamma matrix (first", first.n, " rows&cols) looks like:"))
    print(signif((mypoly$gamma[1:first.n,(1:(first.n*2))]),3))
    print(paste("u vector (first", first.n, "entries) looks like:"))
    print(signif(mypoly$u[1:first.n]))
}

## #' Generic for functions that combine same-type things.
## #' @param obj object
## smartadd <- function(obj, ...) {UseMethod("smartadd")}

##' Adds new polyhedron to an original polyhedron, but only if the inference
##' changes (i.e. Vup and Vlo changes)
##' @param orig.poly Original polyhedron
##' @param new.poly New polyhedron to add, but only if the Vup or Vlo change.
##' @param v contrast vector
##' @param y data vector
smartadd <- function(orig.poly, new.poly, v, y){

    if(all(is.na(orig.poly$gamma[1,])) & nrow(orig.poly$gamma)==1){
        return(new.poly)
    }

    ## Get original Vup and Vlo
    pobj = poly.pval(y = y,
                     G = orig.poly$gamma,
                     u = orig.poly$u,
                     v = v,
                     sigma = sigma)
    Vup = pobj$vup
    Vlo = pobj$vlo


    ## Intersect the original and new polyhedra
    updated.gamma = rbind(orig.poly$gamma, new.poly$gamma)
    updated.u = c(orig.poly$u, new.poly$u)
    updated.poly = structure(list(gamma = updated.gamma, u= updated.u), class = "polyhedra")

    ## Create new Vup and Vlo
    pobj = poly.pval(y = y,
                     G = updated.poly$gamma,
                     u = updated.poly$u,
                     v = v,
                     sigma = sigma)
    Vup.new = pobj$vup
    Vlo.new = pobj$vlo

    ## Compare
    if(Vlo < Vlo.new | Vup > Vup.new ){
        print("Updated polyhedron!")
        return(updated.poly)
    } else {
        return(orig.poly)
    }
}


