## Check uniformity
source("../main/artificial/artif-helpers.R")
n = 10
mn = rep(0,n)
sigma = 1
sigma.add = 0.2
nsim = 500
results = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    y = mn + rnorm(n,0,sigma)
    ## Randomized
    pvs = do_rfl_inference(y=y, max.numSteps=8,
                           consec=2, sigma=sigma, postprocess=TRUE,
                           locs=1:length(y), numIS=100, sigma.add = sigma.add, bits=1000,
                           inference.type="rows")
    ## Nonrandomized inference
    h = binSeg_fixedSteps(y, numSteps = 8)
    consec=2
    ic_obj = get_ic(h$cp, h$y, consec=consec, sigma=sigma, type="bic")
    stoptime = ic_obj$stoptime
    ic.poly = ic_obj$poly
    poly = polyhedra(h, numSteps=stoptime)
    cp = h$cp[1:stoptime]
    cp.sign = h$cp.sign[1:stoptime]
    poly$gamma = rbind(poly$gamma, ic.poly$gamma)
    poly$u = c(poly$u, ic.poly$u)
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    lapply(vlist, function(v){
        poly.pval2(y=y,v=v,poly=poly, sigma=sigma)$pv
    })


    return(pvs)
}, mc.cores=4)

res = results[sapply(results, function(myresult){length(myresult)>1})]

qqunif(unlist(lapply(res, function(myres) myres$pv)))



## Not uniform. What is going wrong? Does fixed inference work? Is there something conceptually wrong?
## Maybe just one step size off?


## See if the inference is correct:
