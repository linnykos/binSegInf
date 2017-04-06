##' Wild binary segmentation, with fixed threshold
##' @param y numeric vector to contain data
##' @param numSteps numeric of number of stepsemas
##' @param numIntervals number of random intervals to sample
##' @param tol tolerance to handle divide by zero instances
##'
##' @return either an environment address which contains output
##'     (\code{env$intervals} and \code{env$tree}), or a list of wild binary
##'     segmentation output.
##' @export
wildBinSeg_fixedSteps <- function(y, numSteps, numIntervals = NULL,
                return.env=FALSE, seed=NULL, verbose=FALSE, intervals = NULL){
    ## Basic checks
    if(numSteps > length(y)-1) stop(paste("You should ask for less than", length(y), "steps!"))
    if(any(duplicated(y))) stop("y must contain all unique values")
    if(is.null(numIntervals) & is.null(intervals)){
        stop("Provide input for generating intervals, or the intervals themselves!")}
    if(!is.null(intervals)) stopifnot(.is_valid_intervals(intervals))
    
    ## Generate the random intervals
    if(!is.null(seed)) set.seed(seed)
    if(is.null(intervals)){
        intervals = generate_intervals(length(y), numIntervals)
    }
    
    ## De-duplicate the intervals
    intervals = .deduplicate_intervals(length(y), intervals)
    starts = intervals$starts
    ends = intervals$ends

    ## Initialize things
    A = T = c()
    S = E = B = Z = M = list()
    Scurr = Ecurr = Bcurr = Zcurr = Mcurr = cplist(2*numSteps)
    Tcurr=Acurr=list()
    
    ## At step 1,
    Scurr = add(Scurr,1,1,1)
    Ecurr = add(Ecurr,1,1,length(y))
    Tcurr[[1]] = c(1,1)
    G = matrix(NA, ncol = length(y), nrow = 2*length(y)*numSteps)
    
    ## At general step
    for(mystep in 2:(numSteps+1)){
        if(mystep==4) browser()

        ## Goal is to get max.m, max.b, max.j, max.k
        .get_max.mbc <- function(Tcurr.without.null){
            
            ## Basic check
            if(any(is.na(Tcurr.without.null))) stop("Input Tcurr should be rid of nulls!")

            ## Get maximizing quantities
            max.m.b.cusums <- lapply(Tcurr.without.null, function(t){
                s = extract(Scurr,t[1],t[2])
                e = extract(Ecurr,t[1],t[2])
                ms = which(.get_which_qualify(s,e,intervals))
                if(length(ms)==0) return()

                ## Get the maximizer (m,b,z)
                max.m = ms[which.max(sapply(ms, function(m) .get_max_b(m,intervals,y,"cusums")))]
                max.cusum = max(sapply(ms, function(m) .get_max_b(m,intervals,y,"cusums")))
                max.z = .get_max_b(max.m, intervals,y,"max.z")
                max.b = .get_max_b(max.m, intervals,y,"max.b")

                return(list(max.m = max.m, max.b = max.b, max.cusum = max.cusum))})
        }

        ## Rid of the null element (the first one)
        Tcurr.without.null <- Tcurr[!sapply(Tcurr, is.null)]
        mbc.list <- .get_max.mbc(Tcurr.without.null)
        if(all(sapply(mbc.list, is.null))) next
        
        ## Extract m,b,z,j,k 
        ind <- which.max(lapply(mbc.list, function(a) if(is.null(a)) FALSE else a$max.cusum))
        jk.max <- Tcurr.without.null[[ind]]
        j.max <-  jk.max[1]
        k.max <-  jk.max[2]
        m.max <- mbc.list[[ind]]$max.m
        b.max <- mbc.list[[ind]]$max.b
        z.max <- sign(mbc.list[[ind]]$max.cusum)

        ## if(verbose) cat("At step", mystep, ", changepoint", b.max, "enters! from (s,e)=",s.max, e.max, fill=TRUE)

        ## Update S and E
        s.max <- extract(Scurr,j.max,k.max)
        e.max <- extract(Ecurr,j.max,k.max)

        ## Change all other *Curr things 
        Bcurr <- add(Bcurr, j.max, k.max, b.max)
        Zcurr <- add(Zcurr, j.max, k.max, z.max)
        Mcurr <- add(Mcurr, j.max, k.max, m.max)

        ## Change active and terminal set
        Acurr[[mystep]] <- c(j.max,k.max)
        Tcurr <- Tcurr.without.null[-ind]
        Tcurr[[mystep]] <- c(j.max + 1, 2*k.max-1)
        Tcurr[[mystep+1]] <- c(j.max + 1, 2*k.max)

        ## Take snapshot
        A[[mystep]] = trim(Acurr)
        T[[mystep]] = trim(Tcurr)
        S[[mystep]] = trim(Scurr)
        E[[mystep]] = trim(Ecurr)
        B[[mystep]] = trim(Bcurr)
        Z[[mystep]] = trim(Zcurr)
        M[[mystep]] = trim(Mcurr)

        ## Prep Scurr and Ecurr for next step.
        Scurr <- add(Scurr, j.max+1, 2*k.max-1, s.max)
        Ecurr <- add(Ecurr, j.max+1, 2*k.max-1, b.max)
        Scurr <- add(Scurr, j.max+1, 2*k.max, b.max+1)
        Ecurr <- add(Ecurr, j.max+1, 2*k.max, e.max)
        
    }

    names(B) = names(M) = names(Z) = names(A) = 
    names(T) = names(S) = names(E) = paste("step", 0:numSteps)

    print(B)
    ## Bundle 
    obj <- structure(list(A=A, T=T, S=S, E=E, B=B, Z=Z, M=M,
                          intervals=intervals, cp=(B[[length(B)]]) [,"val"],
                          cp.sign=(Z[[length(Z)]])[,"val"], numSteps=numSteps, y=y),
                     class="wbsFs")
    return(obj)
}




##' Helper function to trim tree and deduplicate env$signs
##' @param env An environment created as a result of the outer-most run of
##'     \code{.wbs_inner(.., s=1,e=length(y))}
.clean_env <- function(env){
    ## Rid env$tree of the empty element
    env$tree = env$tree[lapply(env$tree,length)>1]

    ## Deduplicate env$signs
    all.m = unique(as.numeric(names(env$signs)))
    unique.first.m.ind  = sapply(all.m, function(my.m){
        min(which(as.numeric(names(env$signs))==my.m))
    })
    env$signs = (env$signs)[unique.first.m.ind]
}

#' is_valid for wbs
#'
#' @param obj wbs object
#'
#' @return TRUE if valid
#' @export
is_valid.wbsFs <- function(obj){
    ## if(!all(names(obj) %in%
    ##         c("tree",
    ##           "y",
    ##           "thresh",
    ##           "signs",
    ##           "intervals",
    ##           "cp",
    ##           "cp.sign",
    ##           "fixedInterval"))) stop("obj must contain certain elements!")
  TRUE
}


## Print wild binary segmentation
print.wbsFs <- function(obj){
    stopifnot(is_valid.wbsFs(obj))
    print("The changepoints of this object are: ")
    print(obj$cp * obj$cp.sign)
}



## #' is_valid for semat; not using this because it destroys the data frame, makes it into a general structure, which screws code.
## #'
## #' @param semat semat object
## #'
## #' @return TRUE if valid
## #' @export
## is_valid.semat <- function(semat){
##     if(!all(colnames(semat) %in% c("m", "b", "maxcusum", "maxhere", "maxhere",
##                                    "passthreshold"))) stop("semat must be a matrix that contains certain elements!") 
##   TRUE
## }

