## First, recovery comparison
onesim_recovery <- function(lev){

    ## Generate data
    n=50
    consec=2
    meanfun=fourjump
    mn = meanfun(lev,n)
    sigma=1
    y = mn + rnorm(n, 0, sigma)
    sigma.add = 0.2
    numSteps=15
    new.noise = rnorm(n,0,sigma.add)
    locs = sapply(1:4, function(ii){ii/5*n + c(-1,0,1)})

    ## Get stopped sbs model
    h = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
    ic_obj = get_ic(h$cp, h$y, consec=consec, sigma=sigma, type="bic")
    stoptime1 = ic_obj$stoptime
    if(is.na(stoptime1)|stoptime1==0){
        cp1 = c()
        correct1 = NA
        error1 = NA
    } else {
        cp1 =h$cp[1:stoptime1]
        correct1 = sum(cp1%in%locs)/length(cp1)
        error1 = (length(cp1) - sum(cp1%in%locs))/length(cp1)
    }

    ## Get stopped wbs model
    g = wildBinSeg_fixedSteps(y + new.noise, numSteps=numSteps, numIntervals=n)
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    stoptime2 = ic_obj$stoptime
    if(is.na(stoptime2)|stoptime2==0) {
        cp1 = c()
        correct2 = NA
        error2 = NA
    } else {
        cp2 = g$cp[1:stoptime2]
        correct2 = sum(cp2%in%locs)/length(cp2)
        error2 = (length(cp2) - sum(cp2%in%locs))/length(cp2)
    }

    ## Get stopped cbs model
    g = circularBinSeg_fixedSteps(y + new.noise, numSteps=numSteps, numIntervals=n)
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    stoptime2 = ic_obj$stoptime
    if(is.na(stoptime2)|stoptime2==0) {
        cp1 = c()
        correct2 = NA
        error2 = NA
    } else {
        cp2 = g$cp[1:stoptime2]
        correct2 = sum(cp2%in%locs)/length(cp2)
        error2 = (length(cp2) - sum(cp2%in%locs))/length(cp2)
    }

    return(data.frame(correct1=correct1,error1=error1,
                      correct2=correct2,error2=error2))
}


## Next, power comparison
onesim_power <- function(lev){

    ## Generate data
    n=50
    consec=2
    meanfun=fourjump
    mn = meanfun(lev,n)
    sigma=1
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    sigma.add = 0.2
    numSteps=15
    numIS=100
    new.noise = rnorm(n,0,sigma.add)
    locs = sapply(1:4, function(ii){ii/5*n + c(-1,0,1)})
    improve.nomass.problem=TRUE
    inference.type="pre-multiply"

    #### 1. Get stopped SBS's power
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
    ic_obj = get_ic(h$cp, h$y, consec=consec, sigma=sigma, type="bic")
    stoptime1 = ic_obj$stoptime
    print('stoptime1 is')
    print(stoptime1)

    ## Get stopped polyhedron
    if(ic_obj$flag=="normal" ){
        ## Get ic-stopped model selection polyhedron
        stopped.gamma = do.call(rbind, g$rows.list[1:(stoptime1+consec)])
        stopped.u = rep(0, nrow(stopped.gamma))
        stopped.poly = polyhedra(obj=stopped.gamma, u=stopped.u)
    } else {
        return(ic_obj$flag)
    }
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=stoptime1+consec)
    poly.fudged = polyhedra(h.fudged)
    ## poly.fudged$gamma = rbind(poly.fudged$gamma, stopped.gamma)
    ## poly.fudged$u = c(poly.fudged$u, stopped.u)

    cp = h$cp[1:stoptime1]
    cp.sign = h$cp.sign[1:stoptime1]

    ## Declutter the changepoints
    do.declutter=TRUE
    if(do.declutter){
        cp = cp[which(cp %in% declutter(cp, how.close = 1, sort=T))]
        cp.sign = cp.sign[which(cp %in% declutter(cp, how.close = 1, sort=T))]
    }

    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)

    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]


    ## Do noise-added inference
    numIS=300
    pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged, stopped.poly = stopped.poly) ## Is this right?
        ## h = binSeg_fixedSteps(y,numSteps=numSteps)
    })
    names(pvs) = (cp*cp.sign)[retain]
    sbs.power = sum(pvs < 0.05/stoptime1)/length(pvs)

    #### 2. Get stopped WBS's power

    ## Fit initial WBS for a generous number of steps
    g = wildBinSeg_fixedSteps(y, numIntervals=n, numSteps=numSteps,
                              inference.type='rows')

    ## Collect the IC information and polyhedron
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    ic_poly = ic_obj$poly

    ## Check for flag
    if(ic_obj$flag=="normal" ){
        ## Get ic-stopped model selection polyhedron
        stopped.gamma = do.call(rbind, g$rows.list[1:(ic_obj$stoptime+consec)])
        stopped.u = rep(0, nrow(stopped.gamma))
        poly = polyhedra(obj=rbind(stopped.gamma, ic_obj$poly$gamma),
                         u=c(stopped.u, ic_obj$gamma$u))
    } else {
        return(ic_obj$flag)
    }
    ## Decide on the stoptime
    stoptime2  = ic_obj$stoptime
    print('stoptime2 is')
    print(stoptime2)

    ## Extract changepoints from stopped model and form contrasts
    cp = g$cp[1:stoptime2]
    cp.sign = g$cp.sign[1:stoptime2]
    better.segment=FALSE
    if(better.segment){
        vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g)
    } else {
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    }
    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Calculate the p-values
    pvs = sapply(vlist, function(v){
        if(randomized){
            cumsum.v = cumsum(v)
            return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                                    sigma=sigma,
                                                    numIS=numIS,
                                                    inference.type=inference.type,
                                                    cumsum.y=cumsum.y,
                                                    cumsum.v=cumsum.v,
                                                    stop.time=stoptime2+consec,
                                                    ic.poly=ic_poly,
                                                    improve.nomass.problem=improve.nomass.problem)
                                                    ))
        } else {
            return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
        }
    })
    names(pvs) = (cp*cp.sign)[retain]
    wbs.power = sum(pvs < 0.05/stoptime2)/length(pvs)
    return(data.frame(wbs.power=wbs.power,sbs.power=sbs.power))
}

