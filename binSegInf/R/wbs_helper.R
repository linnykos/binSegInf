##' Add's a matrix containing information about wbs selection even, as a last
##' element to an existing list of such matrices. The name of this new element
##' is of the form "j,k". This is in lieu of |cplist| which was used to store
##' information for SBS-fixed-threshold events. Also note: more than one
##' maximizer per branch is /fine/!
##' @param mylist a list of matrices
addd = function(mylist,j,k,semat){

    ## Basic checks: list is empty? elements are all correctly formated matrices?
    ## if(length(mylist)==0) stop("You are trying to add to an empty list!")
    ## mylist=NULL
    ## if(length(mylist)==0) mylist= list(1) 
    ## if(!all(sapply(1:length(mylist), function(ii) colnames((mylist)[[ii]])[1]=="m"))) stop("mylist is empty!")

    ## Augment the list with a new element and name it appropriately
    mylist = c(mylist, list(semat))
    names(mylist)[length(mylist)] = paste(j, k, sep = ",")

    ## Return the list
    return(mylist)
}



##' Helper function for a /single/ start and end.
maximize = function(s, e, y, getb=TRUE){

    ## Collect all cusums and signs
    all.bs = (s:(e-1))
    all.cusums = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
    names(all.cusums) = all.bs
    all.cusums = c(rep(NA,s-1), all.cusums, rep(NA,length(y)-s))
    sn = sign(all.cusums) 
    
    ## Obtain cusum-maximizing breakpoint and its sign
    b = which.max(abs(all.cusums))
    z = sn[b]

    ## Return maximizer and sign
    if(getb){ return(b) } else {  return(all.cusums[b]) } 
}



##' Matrix of selection information at this branch
##' @param m indices of intervals that are to be considered (i.e. are contained
##'     in the relevant interval (s:e) at runtime of wbs()
##' @param intervals set of random intervals, drawn between 1 and 60
.make_se_mat = function(m, intervals, y, thresh){
    
    ## Make bare matrix
    mymat = matrix(NA, ncol=5, nrow = length(m), dimnames=NULL)
    mymat = data.frame(mymat)
    colnames(mymat) = c("m", "b", "maxcusum", "maxhere", "passthreshold")
    
    ## Fill in information about selection
    mymat[,1] = m
    mymat[,2] = sapply(intervals$se[m],
                       function(se){ maximize(se[1], se[2], y, TRUE) })
    mymat[,3] = sapply(intervals$se[m],
                       function(se){ maximize(se[1], se[2], y, FALSE) })
    
    ## Mark which guy is the maximizing interval
    mymat[,4] = rep(FALSE,length(m))
    mymat[which.max(abs(mymat[,3])),4] = TRUE
    
    ## Check if that guy's max cusum passed the threshold
    mymat[,5] = rep(FALSE,length(m))
    passed.thresh = (abs(mymat[which.max(abs(mymat[,3])),3]) > thresh)
    mymat[which.max(abs(mymat[,3])),5] = passed.thresh
    
    return(mymat)
}



##' Function to generate random intervals for wild binary segmentation.
##' @param n length of data.
##' @param numIntervals Number of intervals to draw.
##' @param seed seed number for random interval generation; defaults to NULL.
generate_intervals <- function(n, numIntervals, seed=NULL){
    
    ## Basic checks
    if(!is.null(seed)) set.seed(seed)
    
    ## Draw some intervals
    done.drawing = FALSE 
    while(!done.drawing){
        starts = sample(1:n, numIntervals*3, replace=TRUE)
        ends = sample(1:n, numIntervals*3, replace=TRUE)
        reverses = (starts >= ends)
        duplicates = (starts == ends)
        done.drawing = (3*numIntervals-sum(duplicates)> numIntervals)
    }
    
    ## Function to make interval given start and end indices, not necessarily
    ## start<end, with an option to reverse them, or rarely, return NULL when
    ## start == end
    makeInterval = function(start, end, reverse, duplicate){
        if(duplicate) return(NULL)
        if(reverse){
            return(end:start)
        } else {
            return(start:end)
        }
    }
    
    ## Take intervals (identical intervals are eliminated! since they don't play
    ## a role anywhere further along the way, and the max-CUSUM comparisons made
    ## using the de-duplicated set of drawn intervals is still fair/same)
    intervals = Map(makeInterval, starts, ends, reverses, duplicates)
    intervals = intervals[!duplicates]
    intervals = intervals[1:numIntervals]

    ## Startpoints
    starts = sapply(intervals,function(se)se[1])
    ends = sapply(intervals,function(se)se[length(se)])

    ## return
    return(list(starts = starts,
                ends = ends,
                intervals = intervals,
                se = Map(c,starts,ends)))
}



##' Helper function that takes in a tree (list of |semat|'s), and extracts the
##' changepoints
##' @param tree A list of |semat|'s. From the environment |env| you created with
##'     \code{wbs()}, simply use \code{env$tree}.
##' @param returntype One of \code{c("cp","sign")}, for whether to return the
##'     changepoint or the sign of teh changepoints.
.extract_cp_from_tree = function(tree, returntype = c("cp","sign")){
    
    ## Extract changepoints and sign from the tree
    if(returntype == "cp"){ 
        all.cps = sapply(tree,
                         function(semat){
                             passed = which(semat[,"maxhere"] & semat[,"passthreshold"])
                             return(semat[passed,"b"])})
        return(unlist(all.cps))
    } else if (returntype == "sign"){
        all.signs = sapply(tree,
                           function(semat){
                               passed = which(semat[,"maxhere"] & semat[,"passthreshold"])
                               return(sign(semat[passed,"maxcusum"]))})
        return(unlist(all.signs))
    } else {
        stop("|returntype| argument should be one of \"cp\" or \"sign\"")
    }
}




##' Computes the CUSUM (cumulative sum) statistic.
##'
##' Note, we calculate this as the right-to-left difference, by default.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
##' @param right.to.left Whether you want right-to-left difference in the cusum
##'     calculation. Defaults to TRUE.
##' @param contrast.vec If TRUE, then the contrast vector v for cusum=v'y is
##'     returned.
##' @export

cusum <- function(s,b,e,y, right.to.left = TRUE, contrast.vec = FALSE){

    ## Form temporary quantities
    n = e-b+1
    n1 = b+1-s
    n2 = e-b
    if(n==1) stop("n cannot be 1!")
    if(b>=e) stop("b must be strictly smaller than e!")
    if(s>b) stop("b must be larger than or equal to!")

    ## Form contrast
    v = rep(0,length(y))
    v[s:b] = -1/n1 
    v[(b+1):e]  = 1/n2
    v = v * sqrt(1/((1/n1)+(1/n2)))
    if(!right.to.left) v = -v

    ## Return the right thing
    if(contrast.vec) return(v)
    return(sum(v*y))
}

