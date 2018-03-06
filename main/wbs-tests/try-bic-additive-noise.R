source("../main/wbs-tests/sim-helpers.R")


## Figure out how to do BIC x additive noise
nsim=2000
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim)
    n = 20
    consec = 2
    meanfun = onejump
    lev = 0
    mn = meanfun(lev,n)
    sigma = 1
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    sigma.add = 0.2
    numSteps=15
    numIS=100
    new.noise = rnorm(n,0,sigma.add)
    locs = 1:n ## Not conditioning for now.
    improve.nomass.problem = TRUE
    inference.type = "pre-multiply"

    #### 1. Get stopped SBS's power
    h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
    ic_obj = get_ic(h.fudged$cp, y, consec=consec, sigma=sigma, type="bic")
    stoptime = ic_obj$stoptime

    ## Get stopped polyhedron
    if(ic_obj$flag=="normal" ){
        ## Get ic-stopped model selection polyhedron
        poly.fudged = polyhedra(binSeg_fixedSteps(y + new.noise, numSteps=stoptime+consec))
    } else {
        return(ic_obj$flag)
    }

    cp = h.fudged$cp[1:stoptime]
    cp.sign = h.fudged$cp.sign[1:stoptime]
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)

    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add,
                                orig.fudged.poly=poly.fudged, ic.poly=ic_obj$poly)
    })
    names(pvs) = names(vlist)[retain]
    return(pvs)
}, mc.cores=8)

pvs = unlist(results)
pvs = pvs[which(!is.na(pvs))]


## sbs.power = sum(pvs < 0.05/stoptime1)/length(pvs)

    ## Do noise-added inference
    numIS=300
    pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged, stopped.poly = stopped.poly) ## Is this right?
        ## h = binSeg_fixedSteps(y,numSteps=numSteps)
    })
    names(pvs) = (cp*cp.sign)[retain]
    sbs.power = sum(pvs < 0.05/stoptime1)/length(pvs)


## Fixing the IC factor
source("../main/wbs-tests/sim-helpers.R")
## Figure out how to do BIC x additive noise
## Generate data
n=50
consec=2
meanfun=fourjump
lev=5
mn = meanfun(lev,n)
sigma=1
set.seed(0)
y = mn + rnorm(n, 0, sigma)
numSteps=15

#### 1. Get stopped SBS's power
h.fudged = binSeg_fixedSteps(y , numSteps=numSteps)
ic_obj_new = get_ic(h.fudged$cp, y, consec=consec, sigma=sigma, type="bic")

