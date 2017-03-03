##' adds a matrix containing information about wbs selection even, as a last
##' element to an existing list of such matrices. the name of this new element
##' is of the form "j,k". this is in lieu of |cplist| which was used to store
##' information for sbs-fixed-threshold events. also note: more than one
##' maximizer per branch is /fine/!
##' @param mylist a list of matrices
addd = function(mylist,j,k,semat){

    ## Basic checks: list is empty? elements are all correctly formated matrices?
    ## if(length(mylist)==0) stop("you are trying to add to an empty list!")
    ## mylist=null
    ## if(length(mylist)==0) mylist= list(1) 
    ## if(!all(sapply(1:length(mylist), function(ii) colnames((mylist)[[ii]])[1]=="m"))) stop("mylist is empty!")

    ## Augment the list with a new element and name it appropriately
    mylist = c(mylist, list(semat))
    names(mylist)[length(mylist)] = paste(j, k, sep = ",")

    ## return the list
    return(mylist)
}

##' Make an addition to the sign list. Specifically, adds elements
##' @param newsigns the signs to add.
##' @param M the indices in env$intervals
##' @param env environment where signs are embedded in
adddd = function(newsigns,M,env){

    ## Basic checks
    stopifnot(class(newsigns) == "list")

    ## Append the new signs
    Map(function(newsign, m){
        env$signs = c(env$signs,list(newsign))
        names(env$signs)[length(env$signs)] = m
        return()
    },newsigns, M)


    ## No return
    return()
}



##' Helper function for a /single/ start and end.
maximize = function(s, e, y, getb=TRUE){

    ## collect all cusums and signs
    all.bs = (s:(e-1))
    all.cusums = sapply(all.bs, function(b){cusum(s=s, b=b, e=e, y=y)})
    names(all.cusums) = all.bs
    all.cusums = c(rep(NA,s-1), all.cusums, rep(NA,length(y)-s))
    sn = sign(all.cusums) 
    
    ## obtain cusum-maximizing breakpoint and its sign
    b = which.max(abs(all.cusums))
    z = sn[b]

    ## return maximizer and sign
    if(getb){ return(b) } else {  return(all.cusums[b]) } 
   
}



##' Matrix of selection information at this branch
##' @param m indices of intervals that are to be considered (i.e. are contained
##'     in the relevant interval (s:e) at runtime of wbs()
##' @param intervals set of random intervals, drawn between 1 and 60
.make_semat = function(m, intervals, y, thresh){
    
    ## Make bare matrix
    mymat = matrix(NA, ncol=7, nrow = length(m), dimnames=NULL)
    mymat = data.frame(mymat)
    colnames(mymat) = c("m", "s", "b", "e", "maxcusum", "maxhere",
                        "passthreshold")
    
    ## Fill in information about selection
    se = intervals$se[m]
    mymat[,"m"] = m
    mymat[,"s"] = sapply(se, function(vec)vec[1])
    mymat[,"b"] = sapply(intervals$se[m],
                       function(se){
                           maximize(se[1], se[2], y, TRUE) })
    mymat[,"e"] = sapply(se, function(vec)vec[2])
    mymat[,"maxcusum"] = sapply(intervals$se[m],
                       function(se){ maximize(se[1], se[2], y, FALSE) })
    
    ## Mark which guy is the maximizing interval
    mymat[,"maxhere"] = rep(FALSE,length(m))
    mymat[which.max(abs(mymat[,"maxcusum"])),"maxhere"] = TRUE
    
    ## Check if that guy's max cusum passed the threshold
    mymat[,"passthreshold"] = rep(FALSE,length(m))
    max.ind = which.max(abs(mymat[,"maxcusum"]))
    passed.thresh = (abs(mymat[max.ind,"maxcusum"]) > thresh)
    mymat[which.max(abs(mymat[,"maxcusum"])),"passthreshold"] = passed.thresh
    

    return(mymat)
}


##' Matrix of selection information at this branch
##' @param m indices of intervals that are to be considered (i.e. are contained
##'     in the relevant interval (s:e) at runtime of wbs()
##' @param intervals set of random intervals, drawn between 1 and 60
.make_signs = function(m, intervals, y, thresh){
    
    ## Basic checks
    stopifnot(length(m)==1)

    ## Calculate things
    se = intervals$se[[m]]
    s=se[1]; e=se[2];
    all.bs = (s:(e-1))
    all.cusums = sapply(all.bs, function(b){cusum(s=s, b=b, e=e, y=y)})
    names(all.cusums) = all.bs
    sn = sign(all.cusums) 


    ## Make bare matrix
    mymat = matrix(NA, ncol=6, nrow = e-s, dimnames=NULL)
    mymat = data.frame(mymat)
    colnames(mymat) = c("s", "b", "e", "sn", "maxhere", "cusums")
    
    ## Fill in information about selection
    mymat[,"s"] = rep(s, e-s)
    mymat[,"b"] = all.bs
    mymat[,"e"] = rep(e, e-s)
    mymat[,"sn"] = sn
    mymat[, "maxhere"] = FALSE
    mymat[which.max(abs(all.cusums)), "maxhere"] = TRUE
    mymat[,"cusums"] = all.cusums
    
    return(mymat)
}



##' Function to generate random intervals for wild binary segmentation.
##' @param n length of data.
##' @param numIntervals Number of intervals to draw.
##' @param seed seed number for random interval generation; defaults to NULL.
##' @param start.end.list Manual list of starts and ends. Literally an R list with two equal length vectors, each named "start" and "end"
generate_intervals <- function(n, numIntervals, seed=NULL, start.end.list = NULL){
    
    ## Basic checks
    stopifnot(n>1)
    if(!is.null(seed)) set.seed(seed)
    if(!is.null(start.end.list)){
        stopifnot(names(start.end.list) %in% c("start","end"))
        stopifnot(length(start.end.list[["start"]]) == length(start.end.list[["end"]]))
    }

    ## Use inputted start.end.list, or draw intervals
    if(!is.null(start.end.list)){
       starts = start.end.list[[1]]
       ends = start.end.list[[2]]
       reverses = (starts >= ends)
       duplicates = (starts == ends)
       numIntervals = length(starts)
    } else {
        done.drawing = FALSE 
        while(!done.drawing){
            starts = sample(1:n, numIntervals*3, replace=TRUE)
            ends = sample(1:n, numIntervals*3, replace=TRUE)
            reverses = (starts >= ends)
            duplicates = (starts == ends)
            done.drawing = (3*numIntervals-sum(duplicates)> numIntervals)
        }
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
##' @param n length of data. Defaults to \code{length(y)}.
##' @param right.to.left Whether you want right-to-left difference in the cusum
##'     calculation. Defaults to TRUE.
##' @param contrast.vec If TRUE, then the contrast vector v for cusum=v'y is
##'     returned.
##' @param unsigned if TRUE, then returns the contrast vector that makes the
##'     contrast with |y| positive, or the absolute value of the contrast.
##' @export

cusum <- function(s,b,e,n=NULL, y=NULL, right.to.left = TRUE, contrast.vec = FALSE, unsigned = FALSE){

    ## Form temporary quantities
    nn = e-b+1
    n1 = b+1-s
    n2 = e-b
    if(nn==1) stop("n cannot be 1!")
    if(b>=e) stop("b must be strictly smaller than e!")
    if(s>b) stop("b must be larger than or equal to!")
    if(unsigned & is.null(y)) stop("Can't produce an unsigned linear contrast or contrast vector without a y")
    ## if(!is.null(y)) stopifnot(length(y)==n)
    if(is.null(y)){ if(is.null(n)) stop("Must provide one of y or n.")}
    if(is.null(n)) n = length(y)

    ## Form contrast
    v = rep(0,n)
    v[s:b] = -1/n1 
    v[(b+1):e]  = 1/n2
    v = v * sqrt(1/((1/n1)+(1/n2)))
    if(!right.to.left) v = -v
    if(unsigned) v = v * sign(sum(v*y))

    ## Return the right thing
    if(contrast.vec){
        return(v)
    } else if(!is.null(y)) {
        return(sum(v*y))
    } else {
        stop("Check your options for cusum() again!")
    }
}
newcusum = cusum

##' Wrapper for \code{cusum()}, to produce /signed/ cusum (unsigned means that
##' the < contrast , y> is manually made to be positive, and signed means it retains its original
##' sign) contrast given \code{s,b,e} as start, break and end.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
signed_contrast <- function(s,b,e,n=NULL,y){
            cusum(s,b,e,n,y, contrast.vec=TRUE, unsigned=FALSE)
}

##' Wrapper for \code{cusum()}, to produce /unsigned/ cusum contrast given
##' \code{s,b,e} as start, break and end. (Unsigned means that the < contrast ,
##' y> is manually made to be positive, and signed means it retains its original
##' sign).
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
unsigned_contrast <- function(s,b,e,n=NULL,y){
            cusum(s,b,e,n,y, contrast.vec=TRUE, unsigned=TRUE)
}


##' Helper function to /manually/ make contrasts.
##' @param test.bps Breakpoint location to test.
##' @param adj.bps Directly adjacent breakpoint locations.
##' @param sn Sign (direction) of proposed breakpoint at |test.bps|; use +1 or -1.
##' @param n length of data.
##' @examples
##' make_contrast(20,c(1,40),60)
##' @export
make_contrast = function(test.bp, adj.bps, sn, n){

    ## Basic checks
    stopifnot(all(c(test.bp, adj.bps) %in% 1:n))
    stopifnot(min(adj.bps)<=test.bp)
    stopifnot(max(adj.bps)>=test.bp)
    stopifnot(length(sn)==1)
    stopifnot(sn %in% c(-1,1))

    ## Make contrast
    d = rep(0,n)
    d[(min(adj.bps)):(test.bp)] = -1/(test.bp-min(adj.bps)+1)
    d[(test.bp+1):(max(adj.bps))] = +1/(max(adj.bps)-test.bp)
    return(sn*d)
}


##' Helper function for making segment contrasts from a wbs object.
##' @param obj Result from running wbs()
make_all_segment_contrasts <- function(obj){

    ## Basic checks
    stopifnot(is_valid.wildBinSeg(obj))
    if(length(obj$cp)==0) stop("No detected changepoints (by wild binary segmentation)!")

    ## Augment the changepoint set for convenience
    ord = order(obj$cp)
    cp_aug = c(0,obj$cp[ord],n)
    sn_aug = c(NA,obj$cp.sign[ord],NA)
    dlist = list()

    ## Make each contrast
    for(ii in (2:(length(cp_aug)-1))){
        d = rep(0,n)
        ind1 = (cp_aug[ii-1]+1):cp_aug[ii] ## 1 to 3, 4 to 9
        ind2 = (cp_aug[ii]+1):cp_aug[ii+1]
        d[ind1] = -1/length(ind1)
        d[ind2] = 1/length(ind2)
        dlist[[ii-1]] = d * sn_aug[ii]
    }
    names(dlist) = (obj$cp * obj$cp.sign)[ord]
    return(dlist)
}



## ##' Make contrasts from
## ##' @param
## contrast <- function(test.cp, all.cps, all.cps.signs, n){
    
##     ##
##     test.cp = c(4)
##     all.cps = c(4,9,11)
##     all.cps.signs = c(+1,-1,+1)
##     diffs[which.min(diffs)] = +Inf
##     min(diffs)

    
    
##     ik = (obj$pathobjs)$B[k]
##     sk = (obj$pathobjs)$s[k] 
##     breaks = (obj$pathobjs)$B
    
##     if(type == "spike"){
  
##         v = rep(0,length(y))
##         v[ik] = -1
##         v[ik+1] = 1
##         v = sk * v      
        
##     } else if (type == "segment"){
        
##         ## Extract usual segment test endpoints
##         Ks = makesegment(breaks=breaks,k=k,klater=klater,n=length(y))
##         K = Ks$K
##         Kmin = Ks$Kmin
##         Kmax = Ks$Kmax
        
##         ## form vector
##         v = rep(0,length(y))    
##         v[Kmin:K] <- (Kmax - K)/(Kmax - Kmin + 1)
##         v[(K+1):Kmax] <- -(K - Kmin + 1)/(Kmax - Kmin + 1)
##         v <- -sk *v
## }


##' Checking if intervals is correct.
.is_valid_intervals <- function(intervals){
   return(all(names(intervals) %in% c("starts","ends","intervals","se")))
}


##' Deduplicating any intervals.
.deduplicate_intervals <- function(n, intervals){
    ## Basic checks
    stopifnot(.is_valid_intervals(intervals))

    ## Get unique guys, form new intervals
    unique.se = unique(intervals$se)
    unique.start.end.list = list(sapply(unique.se, function(se)se[1]),
                                 sapply(unique.se, function(se)se[2]))
    return(generate_intervals(n=n, start.end.list = unique.start.end.list))
}
