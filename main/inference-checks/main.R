## Synopsis: Conduct inferential checks
source("../main/inference-checks/inference-checks-helpers.R")

sim.settings = list(icstop=TRUE, nulltype=c("zeronull","zeromean")[[1]],
                    type=c("fl-addnoise", "sbs-addnoise-nonmarg", "sbs-addnoise",
                           "wbs", "cbs-addnoise")[[1]])

unifpval_driver <- function(sim.settings){

    ## Fixed settings
    n = numIntervals = 20
    numSteps = 4
    sigma = 1
    sigma.add = 0.2
    meanfun = fourjump
    mean.num.things=20

    ## Conduct simulation; the goal is to just have a vector of p-values in the
    ## end that should be uniform
    lev = (if(sim.settings$nulltype=="zeronull") 0 else 1)

    if(!icstop){
        objs = mclapply(1:nsim, function(isim){
            obj = dosim_compare(type=type, n=n, lev=lev, numIntervals=numIntervals,
                          numIS=10, meanfun=meanfun, numSteps=numSteps,
                          max.numIS=2000, bits=bits,
                          min.num.things=min.num.things,
                          sigma.add=sigma.add, sigma=sigma)
        }, mc.cores=mc.cores)
        if(sim.settings$nulltype=="zeronull"){
            unif.pvs.list = sapply(objs, function(obj)obj$pvs)
            unif.pvs = unlist(unif.pvs.list)
            return(unif.pvs)
        } else {
            ## Filter by v
            v = vlist()

        }
    } else {
        mn = meanfun(lev=lev,n=n)
        y = mn + rnorm(n,0,1)

        ## Obtain inferences.
        how.close = 1
        numIS = 10
        min.num.things = 30
        object = inference_bsFs(y=y, max.numSteps=15, consec=2, sigma=sigma,
                                postprocess=FALSE, locs=1:length(y), numIS= numIS,
                                min.num.things=min.num.things,
                                inference.type="pre-multiply", bits=bits,
                                sigma.add=sigma.add, verbose=verbose,
                                start.time=start.time, how.close=how.close,
                                retain.only.null.cases=FALSE,
                                max.numIS=max.numIS,
                                mn=mn)
    }

    filename <- makesimfilename(sim.settings)

    ## Save the things
    save(unif.pvs, sim.settings, file=file.path(outputdir,filename))
}


## Make a combination of simulation settings
null.types = c("zeronull","zeromean")
inference.types = c("fl-addnoise", "sbs-addnoise-nonmarg", "sbs-addnoise",
                           "wbs", "cbs-addnoise")

## sim.settings = list(icstop=TRUE, null.type=null.types[[1]],
##                     inference.type=infereence.types[[1]])
## sim.settings = list(icstop=FALSE, null.type=null.types[[1]],
##                     inference.type=infereence.types[[1]])
## sim.settings = list(icstop=FALSE, null.type=null.types[[1]],
##                     inference.type=infereence.types[[1]])



