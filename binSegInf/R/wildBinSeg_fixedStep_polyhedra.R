#' Generate polyhedra matrix from wbs output
#' Forms both Gamma matrix and u vector
#'
#' @param obj Output from wbs
#' @param approx If TRUE, then does a Vup/Vlo comparison to see if you should add (chunks) of r
#' @param ... not used now
#'
#' @return An object of class polyhedra
#' @export
polyhedra.wbsFs <- function(obj, approx=FALSE,...){

    ## Basic checks
    stopifnot(is_valid.wbsFs(obj))

    ## Get all polyhedra
    actual.num.steps = (length(obj$B)-1)
    all.steps.polys <- lapply(1:actual.num.steps,
                             function(mystep){ poly_from_snapshot(obj, mystep, approx)})

    ## Combine polyhedron and return
    combined.poly = do.call(combine.polyhedra, all.steps.polys)
    return(combined.poly)
}

##' Takes in object and step, and collects the selection event at that step!
##' This means the characterization of the cusum statistic at (b.max, m.max and
##' s.max) being the largest among all the competing intervals at that step.
##' @param obj \code{wbsFt} object.
##' @param mystep step in the path
##'
##' @return a polyhedra object with the selection event at that step.
poly_from_snapshot <- function(obj, mystep,approx){

    ## Obtain snapshot
    Tcurr <- obj$T[[paste("step",mystep)]]
    Scurr <- obj$S[[paste("step",mystep)]]
    Ecurr <- obj$E[[paste("step",mystep)]]

    ## Extract maximizing things at this step.
    max.m = get_last_row_val(obj$M[[mystep+1]])
    max.b = get_last_row_val(obj$B[[mystep+1]])
    max.z = get_last_row_val(obj$Z[[mystep+1]])

  if(max.m==0){
    jk.max = as.numeric(get_last_row_ind.cplist( (obj$M[[mystep+1]])))
    max.s = extract(Scurr,jk.max[1],jk.max[2])
    max.e = extract(Ecurr,jk.max[1],jk.max[2])
   } else {
    max.s = (obj$intervals$starts)[max.m]
    max.e = (obj$intervals$ends)[max.m]
    }

    ## 1. First, characterize the sign of the max.cusum.contrast
    max.cusum.contrast = unsigned_contrast(max.s, max.b, max.e, y=obj$y)
    newrows1 = rBind(max.cusum.contrast)

    if(!approx){
    ## For each terminal node, characterize the selection event of b.max in m.max.
    newpolylist <- lapply(Tcurr[!sapply(Tcurr,is.null)], function(t){

        ## Get start/end points
        s = extract(Scurr,t[1],t[2])
        e = extract(Ecurr,t[1],t[2])
        ms = which(.get_which_qualify(s,e,obj$intervals))
        if(obj$augment) ms = c(ms,0)
        if(length(ms)==0) return()

        ## 2. Second, Compare /all other/ cusums to that of the grand max
        newrows2 = do.call(rBind, lapply(ms, function(m){
            if(m==0){se = c(s,e)}else{se = obj$intervals$se[[m]]}
            s.to.e = (se[1]:se[2])
            other.bs = s.to.e[-which(s.to.e == se[2])]
            if(m==max.m) other.bs = other.bs[other.bs!=max.b]
            if(length(other.bs)==0) return(rBind(rep(NA,length(obj$y)))[-1,])

            ## Subtract all other contrast from the maximum cusum contrast
            other.cusum.contrasts = do.call(rBind, lapply(other.bs, function(other.b){
                signed_contrast(se[1], other.b, se[2], y=obj$y)}))
            subtracted.contrasts = rBind(sweep(-rBind(other.cusum.contrasts), 2,
                                               max.cusum.contrast, "+" ),
                                         sweep(+rBind(other.cusum.contrasts), 2,
                                               max.cusum.contrast, "+" ))
            if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
            cat("nrow is", nrow(subtracted.contrasts),fill=TRUE)
            return(subtracted.contrasts)
        }))
    }
    if(!approx){

      for(t in Tcurr[!sapply(Tcurr,is.null)]){
        ## Get start/end points
        s = extract(Scurr,t[1],t[2])
        e = extract(Ecurr,t[1],t[2])
        ms = which(.get_which_qualify(s,e,obj$intervals))
        if(obj$augment) ms = c(ms,0)
        if(length(ms)==0) return()

        ## 2. Second, Compare /all other/ cusums to that of the grand max
        newrows2 = do.call(rBind, lapply(ms, function(m){
            if(m==0){se = c(s,e)}else{se = obj$intervals$se[[m]]}
            s.to.e = (se[1]:se[2])
            other.bs = s.to.e[-which(s.to.e == se[2])]
            if(m==max.m) other.bs = other.bs[other.bs!=max.b]
            if(length(other.bs)==0) return(rBind(rep(NA,length(obj$y)))[-1,])

            ## Subtract all other contrast from the maximum cusum contrast
            other.cusum.contrasts = do.call(rBind, lapply(other.bs, function(other.b){
                signed_contrast(se[1], other.b, se[2], y=obj$y)}))
            subtracted.contrasts = rBind(sweep(-rBind(other.cusum.contrasts), 2,
                                               max.cusum.contrast, "+" ),
                                         sweep(+rBind(other.cusum.contrasts), 2,
                                               max.cusum.contrast, "+" ))
            if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
            cat("nrow is", nrow(subtracted.contrasts),fill=TRUE)
            return(subtracted.contrasts)
        }))
        newrows.so.far = rBind(newrows.so.far, newrows)
        ## Continue here -- compare vup and vlo and intersect if things have
        ## changed. Actually, we need to check this across all steps,
        ## /cumulatively/, which means the comparison has to happen /outside/ of
        ## this function, in the main polyhedra.wbsfs() function..
    }

        newrows = rBind(newrows1, newrows2)
        newu = rep(0,nrow(newrows))


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

