#' Wild binary segmentation, with fixed threshold
#' @param y numeric vector to contain data
#' @param numSteps numeric of number of steps
#' @param numInterval number of random intervals to sample
#' @param tol tolerance to handle divide by zero instances
#'
#' @return either an environment address which contains output
#'     (\code{env$intervals} and \code{env$tree}), or a list of wild binary
#'     segmentation output.
#' @export
wbs <- function(y, thresh, numInterval = length(y)*3, return.env=FALSE, seed=NULL){

    ## Basic checks
    if(any(duplicated(y))) stop("y must contain all unique values")
    
    ## Test setting
    ## thresh=2.5
    ## s=1
    ## e=60
    ## j=0; k=1;
    ## set.seed(0)
    ## y = rep(c(0,4),each=30) + rnorm(60,0,1)
    ## numInterval=100
    ## intervals = generate_intervals(s,e,numInterval)

    ## Generate the random intervals
    if(!is.null(seed)) set.seed(seed)

    intervals = generate_intervals(1, n, numInterval)
    starts = intervals$starts
    ends = intervals$ends
    
    ## Create new environment |env|
    env = new.env()
    env$tree = NULL
    env$intervals = intervals ## Carry through the intervals!
    
    ## Run WBS
    .wbs_inner(y, thresh, s, e, 0, 1, verbose, env=env)
    
    ## Rid env$tree of the empty element
    env$tree = env$tree[lapply(env$tree,length)>1]

    wbs.output = list(tree = env$tree,
                      y = y,
                      thresh = thresh,
                      intervals = env$intervals,
                      cp = .extract_cp_from_tree(env$tree, "cp"),
                      cp.sign = .extract_cp_from_tree(env$tree, "sign")
                      )
    
    if(return.env){ return(env) } else{ return(wbs.output) }
}

##' Inner function
.wbs_inner <- function(y, thresh, s, e, j, k, verbose=FALSE, env=NULL){

    cat(j,k,fill=TRUE)

    ## If segment is 2-lengthed, terminate
    if(e-s<=1){
        return() 
        
    ## Otherwise, calculate CUSUMs
    } else {

        ## Collect list of (starts,ends) that qualify at this branch for comparison
        get_which_qualify = function(s,e, intervals){
            return(sapply(intervals$se, function(se){ return(s<=se[1] && se[2]<=e)}))
        }

        ## Extract the qualifying intervals at this branch
        m = which(get_which_qualify(s,e,env$intervals))
        if(length(m)==0) return()

        ## Form a matrix of results and save them
        semat = .make_se_mat(m, env$intervals, y, thresh)
        env$tree = addd(env$tree, j, k, semat)

        ## If threshold is /not/ exceeded, then terminate
        if(all(semat[,"passthreshold"]==F)){
            return()
            
        ## If threshold is exceeded, then recurse.
        } else { 
            ## Extract changepoint
            passed = which(semat[,"maxhere"] & semat[,"passthreshold"])
            b = semat[passed,"b"]
            
            ## Recurse
            .wbs_inner(y, thresh, s,   b,  j+1,  2*k, verbose=verbose, env=env)
            .wbs_inner(y, thresh, b+1, e, j+1, 2*k+1, verbose=verbose, env=env)
        }
    }
}
