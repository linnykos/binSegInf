#' Creates halfspaces for /all/ cusum maximizations that happen, as prescribed
#' in /each/ row of /each/ \code{semat} object in \code{obj$tree}.
##' @param obj Result of running WBS.
max_of_each_obj <- function(obj){
    
    ## Basic check
    stopifnot(is_valid.wildBinSeg(obj))

    ## Apply max to all |semat|'s in |obj|
    polyobjs = lapply(obj$tree, max_of_each_semat, obj) 
    do.call(combine.polyhedra, polyobjs)
}
            

#' Creates halfspaces for the cusum-maximizing /within each/ interval in
#' \code{obj$interval}.
##' @param semat A matrix with results from a node of an obj$tree.
max_of_each_semat <- function(semat, obj){

    ## Form new rows for G and u
    M = semat[,"m"]
    if(length(M)==0)  return(NULL)
    which.are.not.length.two = sapply(M, function(m){
        se=obj$intervals$se[[m]]; return(se[2]-se[1]>=2)})
    M = M[which(which.are.not.length.two)]
    if(length(M)==0)  return(NULL)
    newrows = do.call(rbind, lapply(M, function(m){
        
       ## Get maximum contrast for the m'th interval
        se = obj$intervals$se[[m]]
        b = semat[semat[,"m"]==m,"b"]
        s.to.e = (se[1]:se[2])
        other.bs = s.to.e[-which(s.to.e==b | s.to.e == se[2])]
        max.contrast = unsigned_contrast(s=se[1], b=b, e=se[2], y=obj$y)
        
        ## Subtract all other contrast from the maximum one
        other.cusum.contrasts = do.call(rbind, lapply(other.bs, function(other.b){
            unsigned_contrast(se[1], other.b, se[2], y=obj$y)}))
        subtracted.contrasts = sweep(-rbind(other.cusum.contrasts), 2,
                                    max.contrast, "+" )
        if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
        return(subtracted.contrasts)
    }))
    newu = rep(0,nrow(newrows))

    ## If newrows is empty (no comparisons to be made), then don't do anything
    if(length(as.numeric(newrows))==0){ return(NULL)}
    
    ## Return it as a polyhedron
    return(polyhedra.matrix(obj = newrows, u = newu))
}

#' Creates halfspaces for the main maximize-max-cusums.
##' @param obj Result of running WBS.
max_of_max_obj <- function(obj){
   
    ## Basic check
    stopifnot(is_valid.wildBinSeg(obj))

    ## Apply max to all |semat|'s in |obj|
    ## polyobjs = tryCatch({
    ##     lapply(obj$tree, max_of_max_semat, obj) }, error = function(e){browser()})

    polyobjs = lapply(obj$tree, max_of_max_semat, obj)
    do.call(combine.polyhedra, polyobjs)
}
            
            
##' Creates halfspaces for the main maximize-max-cusums.  TODO: many of the rows
##' overlap; keep track of this and don't add duplicates. Currently, this is not
##' done carefully here; but rather, during aggregation of the polyhedra, the
##' duplicates are eliminated.
##' @param semat A matrix with results from a node of an obj$tree.
max_of_max_semat <- function(semat, obj){

    ## Basic checks
    ## stopifnot(is_valid.semat(semat))
    
    ## Compare the max-cusum of the grand max
    intervals = obj$intervals
    max.ind = semat[,"maxhere"]
    max.m = semat[max.ind,"m"]
    max.s = intervals$se[[max.m]][1]
    max.b = semat[max.ind,"b"]
    max.e = intervals$se[[max.m]][2]

    max.cusum.contrast = newcusum(s=max.s, b=max.b, e=max.e, y=obj$y,
                               contrast.vec=TRUE, unsigned=TRUE)

    ## Form new rows for G and u
    M = semat[which(!max.ind),"m"]
    if(length(M)==0)  return(NULL)
    which.are.not.length.two = sapply(M, function(m){
        se=intervals$se[[m]]; return(se[2]-se[1]>=2)})
    M = M[which(which.are.not.length.two)]
    if(length(M)==0)  return(NULL)
    newrows = do.call(rbind, lapply(M, function(m){
    ## newrows = t(sapply(M, function(m){
        se = intervals$se[[m]]
        b = semat[semat[,"m"]==m,"b"]
        s.to.e = (se[1]:se[2])
        other.bs = s.to.e[-which(s.to.e==b | s.to.e == se[2])]## | s.to.e == se[1] 
        other.cusum.contrasts = do.call(rbind, lapply(other.bs, function(other.b){
            newcusum(se[1], other.b, se[2], n=length(obj$y), y=obj$y,
                  contrast.vec=TRUE, unsigned=TRUE)}))

        subtracted.contrasts = sweep(-rbind(other.cusum.contrasts), 2,
                                     max.cusum.contrast, "+" )
        if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
        return(subtracted.contrasts)
    }))
    newu = rep(0,nrow(newrows))
    
    ## If newrows is empty (no comparisons to be made), then don't do anything
    if(length(as.numeric(newrows))==0){ return(NULL)}
    
    ## Return it as a polyhedron
    return(polyhedra.matrix(obj = newrows, u = newu))
}
    
    
    
##' Creates halfspaces for terminal node handling. Wrapper for \code{check_terminal_semat()}.
##' @param obj Result of running WBS.
check_terminal_obj <- function(obj){

    ## Basic checks
    stopifnot(is_valid.wildBinSeg(obj))
    
    ## Apply terminal check
    polylist = lapply(obj$tree, check_terminal_semat, obj$intervals, obj$thresh, obj$y)
    do.call(combine.polyhedra, polylist)
}




##' Creates halfspaces for terminal node handling.
##' @param semat A matrix with results from a node of an obj$tree.
check_terminal_semat <- function(semat, intervals, thresh, y){

    ## stopifnot(is_valid.semat(semat))
    
    M = semat[,"m"]
    
    ## If this node is a /terminal/ node
    if(all(semat[,"passthreshold"]==F)){
        ## record that /none/ passed the threshold
        newrows = t(sapply(M, function(m){
            se = intervals$se[[m]]
            b = semat[semat[,"m"]==m,"b"]
            other.cusum.contrast = newcusum(s = se[1],
                                         b = b,
                                         e = se[2], n = length(y), y=y,
                                         contrast.vec = TRUE, unsigned=TRUE)
           return(-other.cusum.contrast)
        }))
        newu = rep(-thresh, nrow(newrows))
        ## If this node is /not/ a terminal node
    } else {
        ## record that the grand max passed the threshold
        max.ind = semat[,"maxhere"]
        max.m = semat[max.ind,"m"]
        max.s = intervals$se[[max.m]][1]
        max.b = semat[max.ind,"b"]
        max.e = intervals$se[[max.m]][2]
        max.cusum.contrast = newcusum(max.s, max.b, max.e, n=length(y), y = y,
                                   contrast.vec=TRUE, unsigned=TRUE)
        ## Compare the max-cusum of the grand max to the threshold
        newrows = rbind(max.cusum.contrast)
        newu = thresh
    }        
    return(polyhedra.matrix(obj = newrows, u = newu))
}

##' Creates halfspaces to characterize the signs of cusums
##' @param obj Output from wbs
sign_of_cusums_obj <- function(obj){

    ## Helper to convert a single signmat to new rows
    signmat_into_polyhedra <- function(signmat){
        newrows = do.call(rbind, lapply(1:nrow(signmat), function(irow){
            signmat.row = signmat[irow,]
            newrow = newcusum(s=signmat.row$"s", b=signmat.row$"b", e=signmat.row$"e",
                           n=length(obj$y), contrast.vec=TRUE )
            return(rbind(signmat.row$"sn" * newrow))}))
        return(polyhedra.matrix(obj=newrows, u=rep(0,nrow(newrows))))
    }
    newpoly <- do.call(combine.polyhedra, lapply(obj$signs, signmat_into_polyhedra) )

    ## Return if the polyhedron is correct
    stopifnot(all(newpoly$gamma %*% obj$y > newpoly$u))
    return(newpoly)
}


#' Generate polyhedra matrix from wbs output
#' Forms both Gamma matrix and u vector
#' 
#' @param obj Output from wbs
#' @param ... not use now
#'
#' @return An object of class polyhedra
#' @export
polyhedra.wildBinSeg <- function(obj, ...){
    
    ## Basic checks
    stopifnot(is_valid.wildBinSeg(obj))

    ## Combine polyhedron
    poly = combine.polyhedra(max_of_each_obj(obj),
                             max_of_max_obj(obj),
                             check_terminal_obj(obj),
                             sign_of_cusums_obj(obj))

    return(poly)
}

##' Temporary function to check if it is correct; NOT to be called at
##' runtime, but only at test time or internally.
##' @param poly An object of class polyhedra
##' @param y A data vector with the appropriate length, used when
##'     creating this polyhedron.
check_polyhedra <- function(poly, y){
    ## print(poly$gamma%*%y >= poly$u)
    stopifnot(all(poly$gamma %*% y >= poly$u))
}




## ##' Creates halfspaces for the main maximize-max-cusums. Wrapper for \code{max_of_max_semat()}.
##         intervals = obj$intervals
##         intervals = obj$intervals
##         max.ind = semat[,"maxhere"]
##         max.m = semat[max.ind,"m"]
##         max.s = intervals$se[[max.m]][1]
##         max.b = semat[max.ind,"b"]
##         max.e = intervals$se[[max.m]][2]
##         max.cusum.contrast = cusum(max.s, max.b, max.e, n=length(obj$y), y=obj$y,
##                                    contrast.vec=TRUE, unsigned=TRUE)

        
##         max.ind = semat[,"maxhere"]
##         max.m = semat[max.ind,"m"]
##         max.s = intervals$se[[max.m]][1]
##         max.b = semat[max.ind,"b"]
##         max.e = intervals$se[[max.m]][2]
##         max.cusum.contrast = cusum(max.s, max.b, max.e, n=length(obj$y), y=obj$y,
##                                    contrast.vec=TRUE, unsigned=TRUE)




## aa = lapply(obj$tree, max_of_max_semat, obj)
## semat = obj$tree[[6]]

## for(ii in 1:length(obj$tree)){
## ii=    (length(obj$tree))
## max_of_max_semat(obj$tree[[ii]], obj)
## }
