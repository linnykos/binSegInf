intervals_new <- function(numIntervals, n, comprehensive=FALSE, existing=NULL, distance=0, maxlength=n-1) {

    ## Basic checks
    if(!is.null(existing)){
        assert_that(all(colnames(existing) %in% c("s", "e")),
                   msg = "|existing| must be a 2-column matrix with columns names s and e.")
    }

    ## If |comprehensive| == TRUE, draw /all/ possible intervals.
    if(comprehensive){
        all.se = t(combn(n,2)) ## This is computationally intensive and quite unnecessary.
        x.all = all.se[,1]
        y.all = all.se[,2]
        if(!is.null(existing))stop("Can't do both of: have comprehensive intervals, and exclude |existing| intervals!")
    } else {
        ## numIntervals = 10000
        enough.intervals = FALSE
        x.all = c()
        y.all = c()
        num.valid.intervals = 0
        fac = 5 ## This is just a way to sample a lot initially
        while(!enough.intervals){
            x = sample(n, size=numIntervals*fac, replace=TRUE)
            y = sample(n, size=numIntervals*fac, replace=TRUE)
            too.close = which(unlist(Map(function(a,b){b-a<=distance}, x, y)))
            too.long = which(unlist(Map(function(a,b){b-a>maxlength}, x, y)))
            eliminate = c(too.close, too.long)
            x.all = c(x.all, x[-eliminate])
            y.all = c(y.all, y[-eliminate])
            ## If |existing| matrix is supplied, then exclude these from consideration.
            if(!is.null(existing)){
                to.exclude = unlist(apply(existing,1, function(myrow){
                    which((x.all == myrow["s"]) & (y.all == myrow["e"]))
                }))
                x.all = x.all[-to.exclude]
                y.all = y.all[-to.exclude]
            }
            enough.intervals = (length(x.all) >= numIntervals)
        }
    }

    ## Visualizing the frequency; should be even barplots!
    ## xy = cbind(x.all, y.all)
    ## plot((xy), xlim=c(0,n), ylim=c(0,n))
    ## xylist = lapply(1:nrow(xy), function(irow)xy[irow,])
    ## allpairs = unique(xylist)
    ## freqs = rep(NA,length(allpairs))
    ## for(ipair in 1:length(allpairs)){
    ##     mypair = allpairs[[ipair]]
    ##     freqs[ipair] = sum(sapply(xylist, function(myxy)all(myxy==mypair)))
    ## }
    ## barplot(freqs)

    structure(list(starts=x.all,
                   ends=y.all,
                   cusummat=NA,
                   numIntervals=numIntervals), class="intervals")
}

intervals_old <- function(numIntervals, n, comprehensive=FALSE, existing=NULL, distance=0, maxlength=n-1) {

    ## Basic checks
    if(!is.null(existing)){
        assert_that(all(colnames(existing) %in% c("s", "e")),
                   msg = "|existing| must be a 2-column matrix with columns names s and e.")
    }

    ## Make start-end candidates
    starts = ends = c()
    all.se = t(combn(n,2)) ## This is computationally intensive and quite unnecessary.
    colnames(all.se) = c("s", "e")

    ## Remove the pairs that are too close
    too.close = apply(all.se, 1, function(myrow){
        return((myrow["e"] - myrow["s"]) < distance)
    })
    if(any(too.close)){
        all.se = all.se[-which(too.close),,drop=FALSE]
    }

    ## Remove the pairs that are too long
    too.long = apply(all.se, 1, function(myrow){
        return((myrow["e"] - myrow["s"]) > maxlength)
    })
    if(any(too.long)){
        all.se = all.se[-which(too.long),,drop=FALSE]
    }


    ## If |existing| matrix is supplied, then exclude these from consideration.
    if(!is.null(existing)){
        to.exclude = apply(existing,1, function(myrow){
            which((all.se[,1] == myrow["s"]) & (all.se[,2] == myrow["e"]))
        })
        all.se = all.se[-to.exclude,,drop=FALSE]
    }

    ## If |comprehensive| == TRUE, draw /all/ possible intervals.
    if(comprehensive){
        numIntervals = nrow(all.se)
        random.i = 1:numIntervals
    } else {
        ## Otherwise, draw random intervals, then store the starts and ends.
        random.i = sample(nrow(all.se), numIntervals, replace=FALSE)
    }

    structure(list(starts=all.se[random.i, 1],
                        ends=all.se[random.i, 2],
                        cusummat=NA,
                        numIntervals=numIntervals
                        ), class="intervals")
}



get_rbs_decluttered_cp_model <- function(y=y, max.numSteps=20, numIntervals=length(y),
                              intervals=NULL, consec=2,
                              sigma, postprocess=TRUE, how.close = 5,
                              better.segment=FALSE,
                              locs=1:length(y), numIS=100,
                              inference.type=inference.type,
                              improve.nomass.problem=TRUE, bits=1000,
                              max.numIS=2000,
                              verbose=FALSE, mc.cores=1,
                              min.num.things=30){

    ## Basic checks
    if(!is.null(intervals)){
        if(is.null(intervals)){stop("Provide either |numIntervals| or |intervals|.")}
        numIntervals = intervals$numIntervals
    }

    n = length(y)
    new.noise = rnorm(length(y),0,sigma.add)
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=max.numSteps)
    ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma+sigma.add, type="bic")
    stoptime = ic_obj$stoptime
    if(ic_obj$flag!="normal"){return(ic_obj$flag)}

    if(verbose) cat("stoptime is", stoptime, fill=TRUE)

    ## Check for flag
    if(ic_obj$flag!="normal" ){
        return(NULL)
    }

    ## Extract changepoints from stopped model and declutter
    cp = h.fudged$cp[1:stoptime]
    cp.sign = h.fudged$cp.sign[1:stoptime]
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)

    if(postprocess){
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
    }
    return(cp)
}
