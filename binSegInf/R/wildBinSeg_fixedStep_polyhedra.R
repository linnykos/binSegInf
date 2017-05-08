#' Generate polyhedra matrix from wbs output
#' Forms both Gamma matrix and u vector
#'
#' @param obj Output from wbs
#' @param reduce If TRUE, then does a Vup/Vlo comparison to see if you
#'     should add (chunks) of rows instead of /all/ of them.
#' @param v a contrast vector, if you want to use smart addition of polyhedra.
#' @return An object of class polyhedra
#' @export
polyhedra.wbsFs <- function(obj, v=NULL, reduce=FALSE, sigma=NULL){

    ## Basic checks
    stopifnot(is_valid.wbsFs(obj))
    if(is.null(v) & reduce) stop("Provide v!")
    if(!is.null(v) & is.null(sigma)) stop("Provide sigma!")

    ## Get all polyhedra
    actual.num.steps = (length(obj$B)-1)
    latest.poly = polyhedra(rbind(rep(NA,n)),NA)

    ## Smartly add rows, if the problem size is big
    if(reduce){
        vup = Inf
        vlo = -Inf
        for(mystep in 1:actual.num.steps){
            newpoly = poly_from_snapshot(obj=obj,
                                         mystep=mystep,
                                         reduce=reduce,
                                         v=v,
                                         vup=vup, vlo=vlo,
                                         sigma=sigma)
            vup = newpoly$vup
            vlo = newpoly$vlo
        }
        return(list(v=v,reduce=reduce,obj=obj,vup=vup,vlo=vlo))

    ## Otherwise, just rbind and add all rows!
    } else {
        all.steps.polys <- lapply(1:actual.num.steps,
                                  function(mystep){
            poly_from_snapshot(obj, mystep, reduce)$poly})
        combined.poly = do.call(combine.polyhedra, all.steps.polys)
        return(combined.poly)
    }
}

##' Takes in object and step, and collects the selection event at that step!
##' This means the characterization of the cusum statistic at (b.max, m.max and
##' s.max) being the largest among all the competing intervals at that step.
##' @param obj \code{wbsFt} object.
##' @param mystep step in the path
##' @param reduce If TRUE, then does a Vup/Vlo comparison to see if you should
##'     add (chunks) of the polyhedron or not. This assumes you've chosen the
##'     contrast \code{v} already.
##' @param contrast Contrast vector of interest of inference; only to be used
##'     when \code{reduce==TRUE}.
##'     approximations only work if the checking for change in (Vup, Vlo) is
##'     done /cumulatively/. Defaults to NULL.
##' @param vup vup.
##' @param vlo vlo.
##' @param v v.
##' @param bits bits.
##' @param sigma noise level.
##' @return a polyhedra object with the selection event at that step.
##' @import Matrix
##' @export
poly_from_snapshot <- function(obj, mystep, reduce=FALSE, vup=NULL, vlo=NULL, v=NULL, sigma=NULL, bits=NULL){

    ## Basic checks
    ## if(is.null(latest.poly) & reduce ) stop("Provide the latest polyhedron!")
    ## latest.poly.copy <- latest.poly
    n = length(obj$y)
    if(is.null(vup)) vup = Inf
    if(is.null(vlo)) vlo = -Inf


    ## Obtain snapshot
    Tcurr <- obj$T[[paste("step",mystep)]]
    Scurr <- obj$S[[paste("step",mystep)]]
    Ecurr <- obj$E[[paste("step",mystep)]]

    ## Extract maximizing things at this step.
    max.m = get_last_row_val(obj$M[[mystep+1]])
    max.b = get_last_row_val(obj$B[[mystep+1]])
    max.z = get_last_row_val(obj$Z[[mystep+1]])

    if(max.m==0){
        jk.max = as.numeric(get_last_row_ind.cplist((obj$M[[mystep+1]])))
        max.s = extract(Scurr,jk.max[1],jk.max[2])
        max.e = extract(Ecurr,jk.max[1],jk.max[2])
    } else {
        max.s = (obj$intervals$starts)[max.m]
        max.e = (obj$intervals$ends)[max.m]
    }

    ## 1. First, characterize the sign of the max.cusum.contrast
    max.cusum.contrast = unsigned_contrast(max.s, max.b, max.e, y=obj$y)
    newrows1 = rbind(max.cusum.contrast)

    if(!reduce){
        ## For each terminal node, characterize the selection event of b.max in m.max.

        Tcurr.without.nulls = Tcurr[!sapply(Tcurr,is.null)]

        newpolylist <- lapply(Tcurr.without.nulls, function(t){

            ## Get start/end points
            s = extract(Scurr,t[1],t[2])
            e = extract(Ecurr,t[1],t[2])
            ms = which(.get_which_qualify(s,e,obj$intervals))
            if(obj$augment) ms = c(ms,0)
            if(length(ms)==0) return()

            ## 2. Second, Compare /all other/ cusums to that of the grand max
            newrows2 = matrix(NA, nrow = n*10, ncol = n)
            irow=0

            for(ii in 1:length(ms)){
                m = ms[ii]
                if(m==0){se = c(s,e)}else{se = obj$intervals$se[[m]]}

                s.to.e = (se[1]:se[2])
                other.bs = s.to.e[-which(s.to.e == se[2])]
                if(m==max.m) other.bs = other.bs[other.bs!=max.b]


                if(length(other.bs)!=0){
                    ## Subtract all other contrast from the maximum cusum contrast
                    other.cusum.contrasts = matrix(NA, nrow=length(other.bs), ncol=n)
                    for(jj in 1:length(other.bs)){
                        other.cusum.contrasts[jj,] = signed_contrast(se[1], other.bs[jj], se[2], y=y)
                    }
                    nr = (nrow(other.cusum.contrasts))
                    subtracted.contrasts = matrix(NA, nrow=2*nr, ncol=n)
                    subtracted.contrasts[1:nr,] = sweep(-rbind(other.cusum.contrasts),
                                                       2, max.cusum.contrast, "+" )
                    subtracted.contrasts[(nr+1):(2*nr),] = sweep(+rbind(other.cusum.contrasts), 2,
                                                       max.cusum.contrast, "+" )
                    if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)


                    ## Augment the matrix if needed
                    while(irow + 2*nr > nrow(newrows2)){
                        newrows2 = rbind(newrows2,
                                         matrix(NA, nrow = n*10, ncol = n) )
                    }
                    ## Allocate rows and prep for next step
                    newrows2[irow +(1:(2*nr)), ] = subtracted.contrasts
                    irow = irow + 2*nr

                }
            }
            newrows2 = trim(newrows2)
            newrows = matrix(NA,nrow=nrow(newrows2)+1, ncol = n)
            newrows[1:nrow(newrows2),] = newrows2
            newrows[nrow(newrows2)+1,] = newrows1
            newu = rep(0,nrow(newrows))

            ## If newrows is empty (no comparisons to be made), then don't do anything
            if(length(as.numeric(newrows))==0){ return(NULL)}

            ## Return it as a polyhedron
            return(polyhedra.matrix(obj = newrows, u = newu))
        })
        names(newpolylist) = Tcurr.without.nulls
        vup = NULL
        vlo = NULL
        return(list(poly = do.call(combine.polyhedra, newpolylist),
                    v = v,
                    vup = vup,
                    vlo = vlo,
                    reduce = reduce))
    }

    if(reduce){


        for(t in Tcurr[!sapply(Tcurr,is.null)]){

            ## Get start/end points
            s = extract(Scurr,t[1],t[2])
            e = extract(Ecurr,t[1],t[2])
            ms = which(.get_which_qualify(s,e,obj$intervals))
            if(obj$augment) ms = c(ms,0)
            if(length(ms)==0) next

            ## 2. Second, Compare /all other/ cusums to that of the grand max
            for(m in ms){
                if(m==0){se = c(s,e)}else{se = obj$intervals$se[[m]]}
                s.to.e = (se[1]:se[2])
                other.bs = s.to.e[-which(s.to.e == se[2])]
                if(m==max.m) other.bs = other.bs[other.bs!=max.b]
                if(length(other.bs)==0) next

                ## Subtract all other contrast from the maximum cusum contrast
                other.cusum.contrasts = do.call(rbind, lapply(other.bs, function(other.b){
                                                           signed_contrast(se[1], other.b, se[2], y=obj$y)}))
                subtracted.contrasts = rbind(sweep(-rbind(other.cusum.contrasts), 2,
                                                   max.cusum.contrast, "+" ),
                                             sweep(+rbind(other.cusum.contrasts), 2,
                                                   max.cusum.contrast, "+" ))
                if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
                ## return(subtracted.contrasts)

                ## Add a clump of rows after checking whether Vup & Vlo changes.
                clump.poly = polyhedra(subtracted.contrasts, rep(0,nrow(subtracted.contrasts)))
                vuplo = update_vuplo(poly = clump.poly, v=v, y=obj$y, vup=vup, vlo=vlo, sigma=sigma, bits=bits)
                vlo = vuplo$vlo
                vup = vuplo$vup
            }
        }
        return(list(poly = NULL,
                    v = v,
                    vup = vup,
                    vlo = vlo,
                    reduce = reduce))
    }
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
