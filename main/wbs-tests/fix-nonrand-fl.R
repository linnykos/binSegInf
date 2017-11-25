
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)[8:9]
results.by.lev = list()
mc.cores=8
nsim=3000
nsims=c(seq(from=3000,to=1000,length=5), round(seq(from=600, to=300, length=4) ))[8:9]
visc=visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
n=200

dosim.nonrand.fl <-function(n, lev, numIntervals=n, sigma.add=0.2, numIS=100,
                            meanfun=onejump, visc=NULL, numSteps=1, bits=50){

    type = match.arg(type)
    if(is.null(visc))visc=1:n
    ## numSteps = 1
    sigma = 1
    ## mn = c(rep(0,n/2), rep(lev,n/2))
    mn = meanfun(lev=lev,n=n)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    ## numIS = 100
    inference.type = "pre-multiply"
    improve.nomass.problem = TRUE
    retain=1:n

        ## Get nonrandomized p-value
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f.nonfudged = genlassoinf::dualpathSvd2(y, D=D, maxsteps=4, approx=T)
        Gobj.nonfudged = genlassoinf::getGammat.naive(obj=f.nonfudged, y=y, condition.step=1)
        poly.nonfudged = polyhedra(obj=Gobj.nonfudged$G, u=Gobj.nonfudged$u)
        vlist <- make_all_segment_contrasts(f.nonfudged)
        if(!is.null(visc)){
            retain = which((f.nonfudged$cp %in% visc))
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = (f.nonfudged$cp * f.nonfudged$cp.sign)[retain]

        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })


        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)
        ## Fit binseg on fudged data
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f.fudged = genlassoinf::dualpathSvd2(y+new.noise, D=D, maxsteps=numSteps, approx=T)
        Gobj.fudged = genlassoinf::getGammat.naive(obj=f.fudged, y=y, condition.step=1)
        poly.fudged = polyhedra(obj=Gobj.fudged$G, u=Gobj.fudged$u)

        ## Get randomized p-value
        vlist <- make_all_segment_contrasts(f.fudged)
        if(!is.null(visc)){
            retain = which((f.fudged$cp %in% visc))
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = (f.fudged$cp * f.fudged$cp.sign)[retain]

        pvs.rand = sapply(vlist, function(v){
        pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                sigma.add=sigma.add, orig.fudged.poly= poly.fudged, bits=bits)
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

    if(type=="fl.nonrand"){

        ## Get nonrandomized p-value
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f.nonfudged = genlassoinf::dualpathSvd2(y, D=D, maxsteps=numSteps, approx=T)
        Gobj.nonfudged = genlassoinf::getGammat.naive(obj=f.nonfudged, y=y, condition.step=1)
        poly.nonfudged = polyhedra(obj=Gobj.nonfudged$G, u=Gobj.nonfudged$u)
        vlist <- make_all_segment_contrasts(f.nonfudged)
        if(!is.null(visc)){
            retain = which((f.nonfudged$cp %in% visc))
            if(length(retain)==0){
                return(data.frame(pvs=NA, locs=NA))
            }
            vlist = vlist[retain]
        }
        locs = (f.nonfudged$cp * f.nonfudged$cp.sign)[retain]

        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=y, poly=poly.nonfudged, v=v, sigma=sigma, bits=bits)$pv
        })

        return(data.frame(pvs=pvs,
                          locs=locs))
    }

}
