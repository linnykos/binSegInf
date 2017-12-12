## Figure out how to do BIC x additive noise
    ## Generate data
    n=50
    consec=2
    meanfun=fourjump
lev=1
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

    cp = h.fudged$cp[1:stoptime1]
    cp.sign = h.fudged$cp.sign[1:stoptime1]

    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)

    ## Retain only the changepoints we want results from:
    retain = which((cp %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

v = vlist[[1]]
numIS = 300
pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                        sigma.add=sigma.add,
                        orig.fudged.poly=poly.fudged,
                        stopped.poly=stopped.poly) ## Is this right?

    ## Do noise-added inference
    numIS=300
    pvs = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged, stopped.poly = stopped.poly) ## Is this right?
        ## h = binSeg_fixedSteps(y,numSteps=numSteps)
    })
    names(pvs) = (cp*cp.sign)[retain]
    sbs.power = sum(pvs < 0.05/stoptime1)/length(pvs)
