#' Wild binary segmentation (with fixed steps)
#'
#' y must not have duplicated values. This is to avoid degenerate behavior of
#' binary segmentation. The reported \code{y.fit} is set the lambda where the
#' (k+1)th jump is about to appear.
#'
#' @param y numeric vector to contain data
#' @param numSteps numeric of number of steps
#' @param numInterval number of random intervals to sample
#' @param tol tolerance to handle divide by zero instances
#'
#' @return a flFs object
#' @export
wild_binseg <- function(y, numSteps, numInterval = length(y)*3, tol = 1e-7){
    if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
    if(any(duplicated(y))) stop("y must contain all unique values")
    
    ## Generate the random intervals
    intervals = generate_intervals(1, n, numIntervals)
    starts = intervals$starts
    ends = intervals$ends
    
    ## Run WBS
    n <- length(y)

    ## Create a tree
    sbs <- function(s, e, j, k, X){

        ## Make node (j,k)
        node <- data.tree::Node$new(paste0(j, "-", k))

        ## If segment is 2-lengthed, terminate
        if(e-s<=1){
            ## Not written yet
        ## Otherwise, calculate CUSUMs
        } else {

            ## Collect all cusums and signs
            all.bs = (s:(e-1))
            all.cusums = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
            names(all.cusums) = all.bs
            all.cusums = c(rep(NA,s-1), all.cusums, rep(NA,n-s))
            sn = sign(all.cusums) 
            
            
            ## Obtain cusum-maximizing breakpoint and its sign
            b = which.max(abs(all.cusums))
            z = sn[b]
        }

        ## Check threshold exceedance, then store
        if(abs(all.cusums[b]) < thresh){
            

            b,z,s,e
            env$Blist = add(env$Blist,j+1,k,b)
            env$Zlist = add(env$Zlist,j+1,k,z)
            env$slist = add(env$slist,j,  k,s)
            env$elist = add(env$elist,j,  k,e)
            return(env)
        } else { 
            if(verbose) cat("the biggest cusum was", all.cusums[b], "which passed the threshold:", thresh,  fill=T)

            env$blist = add(env$blist, j+1,k, b)
            env$Blist = add(env$Blist, j+1,k, b)
            env$zlist = add(env$zlist, j+1,k, z)
            env$Zlist = add(env$Zlist, j+1,k, z)
            env$slist = add(env$slist, j  ,k, s)
            env$elist = add(env$elist, j  ,k, e)
        }

        sbs(s,b,j+1,2k)
        sbs(b+1,e,j+1,2k+1)
    
    

    ## Run WBS


    ## structure
    
       
    
  
  y.fit <- .refit_binseg(y, jumps(tree))
    
  structure(list(tree = tree, y.fit = y.fit, numSteps = numSteps), class = "bsFs")
    
    b
    

    for(steps in 1:(numSteps + 1)){
        idx <- .select_nonactive(n, model.mat$Index)
        a.vec <- .compute_fused_numerator(D, idx, y)
        b.vec <- .compute_fused_denominator(D, idx, model.mat[1:(steps-1),,drop = F])
    
        pos.ratio <- a.vec/(1+b.vec)
        pos.ratio[abs(1 + b.vec) < tol] <- 0
        neg.ratio <- a.vec/(-1+b.vec)
        neg.ratio[abs(-1 + b.vec) < tol] <- 0 

    if(max(pos.ratio) > max(neg.ratio)){
      model.mat[steps,] <- c(idx[which.max(pos.ratio)], 1, max(pos.ratio))
    } else {
      model.mat[steps,] <- c(idx[which.max(neg.ratio)], -1, max(neg.ratio))
    }
  }
  
  y.fit <- .refit_flasso(y, model.mat)

  structure(list(model = model.mat[1:numSteps, ], y.fit = y.fit, 
                 numSteps = numSteps), class = "flFs")
}


#' Function to generate random intervals between start (s) and end (e).
#' @param seed seed number for random interval generation; defaults to NULL.
generate_intervals <- function(s, e, numIntervals, seed=NULL){
    if(is.null(seed)) set.seed(seed)
    
    done.drawing = FALSE 
    while(!done.drawing){
        starts = sample(s:e, numIntervals*3, replace=TRUE)
        ends = sample(s:e, numIntervals*3, replace=TRUE)
        reverses = (starts >= ends)
        duplicates = (starts == ends)
        done.drawing = (3*numIntervals-sum(duplicates)> numIntervals)
    }
    
    ## Function to make interval
    makeInterval = function(start, end, reverse, duplicate){
        if(duplicate) return(NULL)
        if(reverse){
            return(start:end)
        } else {
            return(end:start)
        }
    }
    
    ## Take intervals
    intervals = Map(makeInterval, starts, ends, reverses, duplicates)
    intervals = intervals[!duplicates]
    intervals[1:numIntervals]

    ## Startpoints
    starts = sapply(intervals,function(se)se[1])
    ends = sapply(intervals,function(se)se[2])

    ## return
    return(starts = starts,
           ends = ends,
           intervals = intervals)
}


#' is_valid for wbs
#'
#' @param obj wbs object
#'
#' @return TRUE if valid
#' @export
is_valid.wbs <- function(obj){
  ## if(class(obj$tree)[1] != "Node") stop("obj$tree must a Node")
  ## if(!is.numeric(obj$numSteps)) stop("obj$numSteps must be a numeric")
  ## if(length(.enumerate_splits(obj$tree)) != obj$numSteps) 
  ##   stop("obj$tree and obj$numSteps disagree")
  
  TRUE
}
