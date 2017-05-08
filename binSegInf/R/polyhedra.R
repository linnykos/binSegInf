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

    ## ## OLD: Combine separately and return
    ## newgamma2 = do.call(rbind, lapply(polyobjs, function(mypolyobj) mypolyobj$gamma))
    ## newu2 = as.numeric(do.call(c, lapply(polyobjs, function(mypolyobj) mypolyobj$u)))
    ## newpoly = structure(list(gamma = newgamma, u= newu), class = "polyhedra")

    ## NEW: Combine smartly. Handle better!!
    rownums <- sapply(polyobjs, function(a) nrow(a$gamma))
    nonemptyinds = which(rownums!=0)
    starts = cumsum(c(1,rownums[-length(rownums)]))[nonemptyinds]
    ends = cumsum(rownums)[nonemptyinds]
    rowindlist = Map(function(a,b)a:b, starts, ends)
    newgamma = matrix(NA,nrow=sum(rownums), ncol=ncol(polyobjs[[1]]$gamma))
    newu = rep(NA, sum(rownums))
    for(ii in 1:length(nonemptyinds)){
        newgamma[rowindlist[[ii]],] = polyobjs[[nonemptyinds[ii]]]$gamma
        newu[rowindlist[[ii]]] = polyobjs[[nonemptyinds[ii]]]$u
    }

    newpoly = structure(list(gamma = newgamma, u= newu), class = "polyhedra")
    stopifnot(is_valid.polyhedra(newpoly))
    return(newpoly)
}


print <- function(obj, ...) {UseMethod("print")}
## Prints polyhedra
print.polyhedra <- function(mypoly){
    if(nrow(mypoly$gamma)==0 ){ print("Empty polyhedra object!")
    } else if (all(is.na(mypoly$gamma[1,])) & nrow(mypoly$gamma)==1){
        print("Empty polyhedra object!")
    } else {
    first.n = min(10, nrow(mypoly$gamma), ncol(mypoly$gamma)/2)
    print(paste("Gamma matrix (first", first.n, " rows&cols) looks like:"))
    print(signif((mypoly$gamma[1:first.n,(1:(first.n*2))]),3))
    print(paste("u vector (first", first.n, "entries) looks like:"))
    print(signif(mypoly$u[1:first.n]))
    }
}

## #' Generic for functions that combine same-type things.
## #' @param obj object
## smartadd <- function(obj, ...) {UseMethod("smartadd")}

##' Takes in vup and vlo so far, and returns updated list containing vup and vlo
##' only if the inference should change when intersected with the new
##' polyhedron changes (i.e. vup and vlo changes).
##' @param poly Original polyhedron
##' @param v contrast vector
##' @param y data vector
##' @param vup vup so far.
##' @param vlo vlo so far.
##' @param verbose Whether to print things.
##' @param sigma noise level
##' @param bits Number of decimal points for higher precision calculation of
##'     pvalue. (i.e. precision of Rmpfr)
update_vuplo <- function(poly, v, y, vup, vlo, sigma, verbose=FALSE, bits=100){

    ## Get new vup and vlo
    gamma = poly$gamma
    u = poly$u
    poly = structure(list(gamma = gamma, u=u), class = "polyhedra")
    pobj = poly.pval(y = y,
                     G = poly$gamma,
                     u = poly$u,
                     v = v,
                     sigma = sigma,
                     bits = bits)

    vup.new = pmin(vup, pobj$vup)
    vlo.new = pmax(vlo, pobj$vlo)
    pv.new = poly.pval2(vup=vup.new,vlo=vlo.new,y=y,v=v,bits=bits,sigma=sigma)

    ## Return updated vup and vlo
    return(list(vup = vup.new,
                vlo = vlo.new,
                pv = pv.new))
}


