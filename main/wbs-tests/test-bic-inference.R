## Synopsis: Check uniformity of bic inference.


## Fixed step inference
source("../main/artificial/artif-helpers.R")
n = 10
## mn = rep(0,n)
lev = 2
mn = c(rep(0,n/2), rep(lev,n/2))
sigma = 1
sigma.add = 0.2
nsim = 2000
results = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    y = mn + rnorm(n,0,sigma)
    ## ## Randomized
    ## pvs = do_rfl_inference(y=y, max.numSteps=8,
    ##                        consec=2, sigma=sigma, postprocess=TRUE,
    ##                        locs=1:length(y), numIS=100, sigma.add = sigma.add, bits=1000,
    ##                        inference.type="rows")

    ## Nonrandomized inference
    h = binSeg_fixedSteps(y, numSteps = 8)
    consec = 2
    ic_obj = get_ic(h$cp, h$y, consec=consec, sigma=sigma, type="bic")
    if(ic_obj$flag!="normal") return(NULL)
    stoptime = ic_obj$stoptime
    ic.poly = ic_obj$poly
    poly = polyhedra(h, numSteps=stoptime+consec)
    cp = h$cp[1:stoptime]
    cp.sign = h$cp.sign[1:stoptime]
    poly$gamma = rbind(poly$gamma, ic.poly$gamma)
    poly$u = c(poly$u, ic.poly$u)
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    pvs = lapply(vlist, function(v){
        poly.pval2(y=y,v=v,poly=poly, sigma=sigma)$pv
    })

    truths = lapply(vlist, function(v){
        return(v%*%mn!=0)
    })
    return(pvs)
}, mc.cores=3)

## Nonrandomized fixed inference doesn't seem to have uniform p-values.. Why?
res = results[sapply(results, function(myresult){length(myresult)>1})]
qqunif(unlist(res))
(unlist(lapply(res, function(myres) myres$pv)))

## Slightly /anticonservative/ for fixed p-values..


## Now trying /randomized/ inference

n = 10
## mn = rep(0,n)
lev=2
mn = c(rep(0,n/2), rep(lev,n/2))
sigma = 1
sigma.add = 0.2
nsim = 2000
max.numSteps=8
locs=1:n
inference.type = "pre-multiply"
results = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    y = mn + rnorm(n,0,sigma)
    cumsum.y = cumsum(y)

    ## Fit model and get IC information
    new.noise = rnorm(length(y),0,sigma.add)
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=max.numSteps)
    ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma+sigma.add, type="bic")
    stoptime = ic_obj$stoptime
    if(ic_obj$flag!="normal"){return(ic_obj$flag)}

    ## Stopped fudged model (not used because we can't store the polyhedra)
    ## h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=stoptime+consec)
    ## poly.fudged = polyhedra(h.fudged)

    ## Collect stopped model and postprocess
    cp = h.fudged$cp[1:stoptime]
    cp.sign = h.fudged$cp.sign[1:stoptime]
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y))
    ## if(postprocess){
    ##     cpobj = declutter(cp=cp, cp.sign)$cp.sign
    ##     cp = cpobj$cp
    ##     cp.sign = cpobj$cp.sign
    ## }

    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Do noise-added inference
    pvs = sapply(vlist, function(v){
        cumsum.v = cumsum(v)
        pv = randomize_addnoise(y= y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.obj = h.fudged,
                                numSteps = stoptime+consec,
                                ic.poly = ic_obj$poly, bits=bits,
                                inference.type=inference.type)
        return(pv)
    })
    names(pvs) = names(vlist)

    ## Get truths
    truths = lapply(vlist, function(v){
        return(v%*%mn!=0)
    })
    names(truths) = names(vlist)

    return(pvs)
}, mc.cores=4)

## Nonrandomized fixed inference doesn't seem to have uniform p-values.. Why?
res = results[sapply(results, function(myresult){length(myresult)>1})]
qqunif(unlist(res))
(unlist(lapply(res, function(myres) myres$pv)))
