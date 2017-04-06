#' Generate polyhedra matrix from wbs output
#' Forms both Gamma matrix and u vector
#' 
#' @param obj Output from wbs
#' @param ... not used now
#'
#' @return An object of class polyhedra
#' @export
polyhedra.wbsFs <- function(obj, ...){
    
    ## Basic checks
    stopifnot(is_valid.wbsFs(obj))

    ## Get all polyhedra
    all.step.polys <- lapply(1:maxsteps, function(mystep){ poly_from_snapshot(obj, mystep)})

    ## Combine polyhedron and return
    return(combine.polyhedra(all.steps.polys))
}

##' Takes in object and step, and collects the selection event at that step!
##' This means the characterization of the cusum statistic at (b.max, m.max and
##' s.max) being the largest among all the competing intervals at that step.
##' @param obj \code{wbsFt} object.
##' @param mystep step in the path
##' 
##' @return a polyhedra object with the selection event at that step.
poly_from_snapshot <- function(obj, mystep){
        ## Obtain snapshot
        Tcurr <- obj$T[[mystep]]
        Scurr <- obj$S[[mystep]]
        Ecurr <- obj$E[[mystep]]
    
        ## For each terminal node, characterize the selection event of b.max in m.max.
        newpolylist <- lapply(Tcurr, function(t){
            
            ## Get start/end points
            s = extract(Scurr,t[1],t[2])
            e = extract(Ecurr,t[1],t[2])
            ms = which(.get_which_qualify(s,e,intervals))
            if(length(ms)==0) return()
            
            ## Extract the snapshot.
            max.m = obj$M[[mystep]]
            max.b = obj$B[[mystep]]
            max.z = obj$Z[[mystep]]
            
            ## 1. First, characterize the sign of the max.cusum.contrast
            max.cusum.contrast = unsigned_contrast(max.s, max.b, max.e, y=obj$y)
            newrows1 = rbind(max.cusum.contrast)
            newu1 = rep(0,nrow(newrows1))
            
            ## 2. Second, Compare the max-cusum of the grand max
            newrows2 = do.call(rbind, lapply(ms, function(m){
                se = intervals$se[[m]]
                b = semat[semat[,"m"]==m,"b"]
                s.to.e = (se[1]:se[2])
                other.bs = s.to.e[-which(s.to.e == se[2])] 
                if(m==max.m) other.bs = other.bs[other.bs!=max.b]
                if(length(other.bs)==0) return(rbind(rep(NA,length(obj$y)))[-1,])
                
                ## Subtract all other contrast from the maximum cusum contrast
                other.cusum.contrasts = do.call(rbind, lapply(other.bs, function(other.b){
                    signed_contrast(se[1], other.b, se[2], y=obj$y)}))
                subtracted.contrasts = rbind(sweep(-rbind(other.cusum.contrasts), 2,
                                                   max.cusum.contrast, "+" ),
                                             sweep(+rbind(other.cusum.contrasts), 2,
                                                   max.cusum.contrast, "+" ))
                if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
                return(subtracted.contrasts)
            }))
            newrows = rbind(newrows1, newrows2)
            newu = c(newu1, newu2)
            
            ## If newrows is empty (no comparisons to be made), then don't do anything
            if(length(as.numeric(newrows))==0){ return(NULL)}
            
            ## Return it as a polyhedron
            return(polyhedra.matrix(obj = newrows, u = newu))
        })
    
        return(do.call(combine.polyhedra, newpolylist))
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

