##' Fixed step WBS function. There is a ``mimic'' option that helps you create
##' new polyhedron that says that the same |numSteps| winners win, but among a
##' new set of intervals.
##' @param intervals Optionally, add an |intervals| class object that contains
##'     intervals.
##' @return List of useful algorithm results, including the |gamma| matrix and |u|.
wildBinSeg_fixedSteps <- function(y, numSteps, numIntervals=NULL,
                                  intervals=NULL, mimic=FALSE, wbs.obj=NULL,
                                  comprehensive=FALSE, cumsum.y=NULL, cumsum.v=NULL,
                                  inference.type =c("rows","pre-multiply")){

    inference.type = match.arg(inference.type)

    ## Basic checks
    n = length(y)
    assert_that(!is.null(numIntervals) | !is.null(intervals),
                msg="Must provide either |numIntervals| or |intervals|.")
    if(mimic & is.null(wbs.obj))stop("When mimic=TRUE, provide a wbs object.")
    if(!mimic & !is.null(wbs.obj))stop("Even though mimic==FALSE, you're providing a wbs object.")

    ## Calculate all cusums for all intervals
    if(is.null(intervals)){
        intervals = intervals(numIntervals=numIntervals, n= n,
                              comprehensive=comprehensive)
    }
    intervals = addcusum(intervals=intervals, y=y)

    ## Make empty array to store selection events
    results = matrix(NA, nrow=numSteps, ncol = 6)
    colnames(results) = c("max.s", "max.b", "max.e", "max.sign", "scope.s", "scope.e")

    ## Initially, everything qualifies
    qual.inds = 1:nrow(intervals$cusummat)
    scope = rbind(c(1,n))
    colnames(scope)  = c("s", "e")
    new.rows = new.info = list()

    ## Iterate
    istep = 1
    time.to.stop = FALSE
    while(istep <= numSteps & !time.to.stop){

        ## Do maximization in the respective intervals.
        if(mimic){
            max.info = .get_max_info(wbs.obj, istep)
        } else {
            max.info = maximize.intervals(intervals=intervals, qual.inds=qual.inds)
        }
        if(inference.type=="rows"){
        ## Collect halfspaces for the model selection event
        new.rows[[istep]] = form_rows(intervals=intervals,
                                      max.s=max.info$max.s,
                                      max.b=max.info$max.b,
                                      max.e=max.info$max.e,
                                      max.sign=max.info$max.sign,
                                      qual.inds=qual.inds,
                                      n=n,
                                      y=(if(!mimic) y else NULL))
        } else {

        ## New routine for new.rows, that collects Gy and Gv instead of actual rows
        new.info[[istep]] = form_info(intervals=intervals,
                                      max.s=max.info$max.s,
                                      max.b=max.info$max.b,
                                      max.e=max.info$max.e,
                                      max.sign=max.info$max.sign,
                                      qual.inds=qual.inds,
                                      cumsum.y=cumsum.y,
                                      cumsum.v=cumsum.v)
        }

        ## Calculate the maximizing scope
        max.scope.irow = which(apply(scope, 1,function(myscope){
            (myscope["s"]<= max.info$max.b) & (max.info$max.b <= myscope["e"])
        }))
        max.scope = scope[max.scope.irow,]

        ## Update qualifying indices in cusummat
        old.qual.inds = restrict(intervals=intervals,  s=max.scope["s"], e=max.scope["e"])
        new.qual.inds1 = restrict(intervals=intervals, s=max.scope["s"], e=max.info$max.b)
        new.qual.inds2 = restrict(intervals=intervals, s=max.info$max.b+1, e=max.scope["e"])
        qual.inds = qual.inds[!(qual.inds %in% old.qual.inds)]
        qual.inds = c(qual.inds, new.qual.inds1, new.qual.inds2)

        ## Update scopes
        scope = scope[-max.scope.irow,,drop=FALSE]
        scope = rbind(scope,
                      c(max.scope["s"], max.info$max.b),
                      c(max.info$max.b + 1, max.scope["e"]))
        too.short = which(scope[,"e"] - scope[,"s"] < 1)
        if(length(too.short)>=1){
            scope = scope[-too.short,,drop=FALSE]
        }

        ## Record some results
        results[istep,] <- c(max.info$max.s, max.info$max.b, max.info$max.e, max.info$max.sign,
                             max.scope["s"], max.scope["e"])

        ## Break if done.
        time.to.stop = (nrow(scope)==0 | length(qual.inds)==0)
        istep = istep+1
    }
    stoptime = istep-1
    if(stoptime != numSteps) print("Stopped early!")
    results = results[1:stoptime,,drop=FALSE]

    ## Aggregate things for inference
    if(inference.type=="rows"){
        gamma=do.call(rbind, new.rows)
        names(new.rows) = paste0("step-",1:stoptime)
        Gy=NULL
        Gv=NULL
        u = rep(0,nrow(gamma))
    } else {
        Gy=do.call(c,lapply(new.info, function(a)a[["Gy"]]))
        Gv=do.call(c,lapply(new.info, function(a)a[["Gv"]]))
        gamma=NULL
        u=rep(0, length(Gy))
    }

    return(structure(list(results=results, gamma=gamma, u=u,
                          Gy=Gy, Gv=Gv,
                          rows.list = new.rows,
                          info.list = new.info,
                          cp=results[,"max.b"], cp.sign=results[,"max.sign"],
                          numSteps=numSteps,
                          mimic=mimic,
                          intervals=intervals,
                          numIntervals=numIntervals,
                          y=y), class="wbsFs"))
}

##' Print function for convenience, of |wbs| class object.
print.wbsFs <- function(obj){
    if(obj$mimic) cat("Mimicked object!", fill=TRUE)
    cat("Detected changepoints using WBS with", obj$numSteps, "steps is", obj$cp * obj$cp.sign, fill=TRUE)
}

##' Checks if obj is a valid |wbs| class object
is_valid.wbsFs <- function(obj){
    return(all(names(obj) %in% c("results", "gamma", "u", "cp", "cp.sign", "y",
                                 "numSteps", "mimic", "rows.list","intervals",
                                 "numIntervals", "Gy", "Gv", "info.list")))
}

.get_max_info <- function(wbs.obj,...){ UseMethod(".get_max_info")}
##' Helper to /manually/ obtain maximizing information from a |wbs| class object, at i'th step
.get_max_info.wbsFs <- function(wbs.obj, istep){
    assert_that(is_valid.wbsFs(wbs.obj))
    myrow  = wbs.obj$results[istep,]
    mylist = lapply(myrow, function(a) a)
    mylist = mylist[names(mylist) %in% c("max.s", "max.b", "max.e", "max.sign")]
    mylist[["max.i"]] = NA
    return(mylist)
}




#########################
## Old fixed-step code ##
#########################



## ##' Wild binary segmentation, with fixed threshold
## ##' @param y numeric vector to contain data
## ##' @param numSteps numeric of number of stepsemas
## ##' @param numIntervals number of random intervals to sample
## ##' @param tol tolerance to handle divide by zero instances
## ##' @param intervals Intervals to manually use, produced from
## ##'   \code{generate_intervals()}. Defaults to NULL.
## ##' @param augment TRUE if the qualifying set of intervals in wild binary
## ##'   segmentation should be augmented by the large interval \code{c(s,e)}.
## ##'   Defaults to FALSE.
## ##'
## ##' @return either an environment address which contains output
## ##'     (\code{env$intervals} and \code{env$tree}), or a list of wild binary
## ##'     segmentation output
## ##' @export
## wildBinSeg_fixedSteps <- function(y, numSteps, numIntervals = NULL,
##                                       return.env=FALSE, seed=NULL, verbose=FALSE,
##                                       intervals = NULL, augment=TRUE){
##     ## Basic checks
##     if(numSteps > length(y)-1) stop(paste("You should ask for less than", length(y), "steps!"))
##     if(round(numSteps) != numSteps) stop(paste("You should provide an integer value for numSteps!"))
##     ## if(any(duplicated(y))) stop("y must contain all unique values")
##     if(is.null(numIntervals) & is.null(intervals)){
##         stop("Provide input for generating intervals, or the intervals themselves!")}
##     if(!is.null(intervals)) stopifnot(.is_valid_intervals(intervals))

##     ## Generate the random intervals
##     if(!is.null(seed)) set.seed(seed)
##     if(is.null(intervals)){
##         intervals = generate_intervals(length(y), numIntervals)
##     }
##     intervals = .deduplicate_intervals(length(y), intervals)

##     ## Pre-calculate things
##     cumsums <- cumsum(y)
##     info <- Map(function(s,e){
##         get_morethan_cusums2(s,e,cumsums)},
##         intervals$starts, intervals$ends)

##     ## Initialize things
##     A = T = c()
##     S = E = B = Z = M = list()
##     Scurr = Ecurr = Bcurr = Zcurr = Mcurr = cplist(2*numSteps)
##     Tcurr = Acurr = list()

##     ## At step 1,
##     Scurr = add(Scurr,1,1,1)
##     Ecurr = add(Ecurr,1,1,length(y))
##     Tcurr[[1]] = c(1,1)
##     G = matrix(NA, ncol = length(y), nrow = 2*length(y)*numSteps)

##     ## At general step
##     for(mystep in 2:(numSteps+1)){

##         ## Goal is to get max.m, max.b, max.j, max.k
##         .get_max.mbc <- function(Tcurr){

##           ## Get maximizing quantities
##           max.m.b.cusums <- lapply(Tcurr, function(tt){
##               if(is.null(tt)) return()
##               s = extract(Scurr,tt[1],tt[2])
##               e = extract(Ecurr,tt[1],tt[2])
##               ms = which(.get_which_qualify(s,e,intervals))
##               if(augment) ms = c(ms,0)
##               if(length(ms)==0) return()

##               ## Augmentation is painful.. programmatically
##               info.aug = c(info, list(get_morethan_cusums2(s,e,cumsums)))
##               ms.aug = ms
##               ms.aug[which(ms.aug==0)] = length(intervals$se) + 1

##               ## Get the maximizer (m,b,z)
##               maxcusums = sapply(info.aug[ms.aug], '[[', "max.cusum")
##               max.ind.of.maxcusums = which.max(abs(maxcusums))
##               max.m = ms[max.ind.of.maxcusums]
##               max.cusum = maxcusums[max.ind.of.maxcusums]
##               max.z = sapply(info.aug[ms.aug], '[[', "max.z")[max.ind.of.maxcusums]
##               max.b = sapply(info.aug[ms.aug], '[[', "max.b")[max.ind.of.maxcusums]
##               ## if(tt[1]==2 & tt[2]==1 ) print(c(s,e,max.b,max.z, max.cusum))

##               return(list(max.m = max.m, max.b = max.b, max.cusum = max.cusum, max.z=max.z))})
##         }
##       mbc.list <- .get_max.mbc(Tcurr)
##       if(all(sapply(mbc.list, is.null))){
##         if(verbose){
##           print(paste("There were no qualifying intervals during step ", mystep-1));
##         }
##         next
##       }

##       ## Extract m,b,z,j,k
##       ind <- which.max(lapply(mbc.list, function(a) if(is.null(a)){ -Inf} else{ abs(a$max.cusum)}))
##         ## ind <- which.max(lapply(mbc.list, function(a) if(is.null(a)) FALSE else a$max.cusum))


##       ## ind <- which.max(lapply(mbc.list, function(a) if(is.null(a)){ -Inf} else{ abs(a$max.cusum)}))
##       jk.max <- Tcurr[[ind]]
##       j.max <-  jk.max[1]
##       k.max <-  jk.max[2]
##       m.max <- mbc.list[[ind]]$max.m
##       b.max <- mbc.list[[ind]]$max.b
##       z.max <- mbc.list[[ind]]$max.z  ## sign(mbc.list[[ind]]$max.cusum)
##       ## if(j.max==2 & k.max==1)  browser()

##       ## Update S and E
##       s.max <- extract(Scurr,j.max,k.max)
##       e.max <- extract(Ecurr,j.max,k.max)

##       ## if(verbose) cat("at step", mystep, ", changepoint", b.max,
##       ##                 "enters! from (s,e)=",intervals$se[[m.max]], fill=TRUE)
##       ## if(verbose) cat("at step", mystep, ", changepoint", b.max, "enters! from (s,e)=",s.max,e.max, fill=true)

##       ## Change all other *Curr things
##       Bcurr <- add(Bcurr, j.max, k.max, b.max)
##       Zcurr <- add(Zcurr, j.max, k.max, z.max)
##       Mcurr <- add(Mcurr, j.max, k.max, m.max)

##       ## Take snapshot
##       A[[mystep]] = trim(Acurr)
##       T[[mystep]] = trim(Tcurr)
##       S[[mystep]] = df_to_cplist(trim(Scurr))
##       E[[mystep]] = df_to_cplist(trim(Ecurr))
##       B[[mystep]] = df_to_cplist(trim(Bcurr))
##       Z[[mystep]] = df_to_cplist(trim(Zcurr))
##       M[[mystep]] = df_to_cplist(trim(Mcurr))

##       ## Change active and terminal set
##       Acurr[[mystep]] <- c(j.max,k.max)
##       Tcurr <- Tcurr[-ind]
##       Tcurr[[mystep]] <- c(j.max + 1, 2*k.max-1)
##       Tcurr[[mystep+1]] <- c(j.max + 1, 2*k.max)


##       ## Prep Scurr and Ecurr for next step.
##       Scurr <- add(Scurr, j.max+1, 2*k.max-1, s.max)
##       Ecurr <- add(Ecurr, j.max+1, 2*k.max-1, b.max)
##       Scurr <- add(Scurr, j.max+1, 2*k.max, b.max+1)
##       Ecurr <- add(Ecurr, j.max+1, 2*k.max, e.max)

##       ## Prune Tcurr of all one-length segments
##       Tcurr = prune_of_1_length_segments(Tcurr,Scurr,Ecurr)
##     }

##     if(numSteps==1) A =list("step 0" = NULL, "step 1" = NULL);
##     names(B) = names(M) = names(Z) =
##     names(T) = names(S) = names(E) = paste("step", 0:(length(B)-1))

##     ## Temporary addition of winning intervals, for one-step WBS.
##     if(numSteps==1){
##         winning.interval = c(s.max, e.max)
##     } else {
##         winning.interval = NULL
##     }

##     ## Bundle
##     obj <- structure(list(A=A, T=T, S= S, E=E, B=B, Z=Z, M=M,
##                           intervals=intervals,
##                           winning.interval = winning.interval,
##                           cp=((B[[length(B)]])$mat) [,"val"],
##                           cp.sign=((Z[[length(Z)]])$mat)[,"val"],
##                           numSteps=numSteps, y=y, augment=augment),
##                      class="wbsFs")
##     return(obj)
## }



## ## Helper function to trim tree and deduplicate env$signs, for WBS-ft and WBS-fs
## ## @param env An environment created as a result of the outer-most run of
## ##     \code{.wbs_inner(.., s=1,e=length(y))}
## .clean_env <- function(env){
##     ## Rid env$tree of the empty element
##     env$tree = env$tree[lapply(env$tree,length)>1]

##     ## Deduplicate env$signs
##     all.m = unique(as.numeric(names(env$signs)))
##     unique.first.m.ind  = sapply(all.m, function(my.m){
##         min(which(as.numeric(names(env$signs))==my.m))
##     })
##     env$signs = (env$signs)[unique.first.m.ind]
## }


## #' is_valid for wbs
## #'
## #' @param obj wbs object
## #'
## #' @return TRUE if valid
## #' @export
## is_valid.wbsFs <- function(obj){
##     ## if(!all(names(obj) %in%
##     ##         c("tree",
##     ##           "y",
##     ##           "thresh",
##     ##           "signs",
##     ##           "intervals",
##     ##           "cp",
##     ##           "cp.sign",
##     ##           "fixedInterval"))) stop("obj must contain certain elements!")
##   TRUE
## }


## ## Print wild binary segmentation
## print.wbsFs <- function(obj){
##     stopifnot(is_valid.wbsFs(obj))
##     print("The changepoints of this object are: ")
##     print(obj$cp * obj$cp.sign)
## }



## ## #' is_valid for semat; not using this because it destroys the data frame, makes it into a general structure, which screws code.
## ## #'
## ## #' @param semat semat object
## ## #'
## ## #' @return TRUE if valid
## ## #' @export
## ## is_valid.semat <- function(semat){
## ##     if(!all(colnames(semat) %in% c("m", "b", "maxcusum", "maxhere", "maxhere",
## ##                                    "passthreshold"))) stop("semat must be a matrix that contains certain elements!")
## ##   TRUE
## ## }
