## Synopsis: Checking the uniformity of three methods (1) noise-added binseg,
## (2) fixed binseg, and (3) WBS (randomized) -- in terms of two things (a)
## delta=0, one and two steps (b) delta>0 for several steps testing contrasts
## whose means are zero.


## Simulation settings
n = 20
levs = c(0, 1)
results.by.lev = list()
mc.cores = 6
nsim = 1000
visc = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
numIS = 100
meanfun = fourjump
bits=1000

## (a) For flat signal, collect three things
mylev = 0

## (1a) All three methods
onething <- function(lev){
    types = c("sbs.rand", "sbs.nonrand", "wbs.rand")
    pvs.by.method = sapply(1:length(types), function(itype){
        ## itype = 1
        mytype = types[[itype]]
        visc = 1:n
        results = list()
        for(numsteps in c(1,2)){
            myresult = mclapply(1:nsim, function(isim){
                printprogress(isim,nsim)
                dosim_compare(type="sbs.nonrand", n=n, lev=lev, numIS=numIS, meanfun=meanfun,
                              visc=visc, numSteps=numSteps, bits=bits,
                              max.numIS=max.numIS)
            })
            results[[paste0("Step ", numsteps)]] = myresult
        }

        ## For flat
        pvs.list = sapply(results, function(myresult){
            pvs = sapply(myresult, function(a) a$pvs)
            return(pvs)
        })

        ## For nonflat things TODO: make this so it only collects nonnull
        ## contrasts! this may take a little more coding time
        pvs.list = sapply(results, function(myresult){
            pvs = sapply(myresult, function(a) a$pvs)
            return(pvs)
        })


        names(pvs.list) = names(results)
    })
    names(pvs.by.method) = types
    return(pvs.by.method)
}

flat.pvs.list = onething(0)

## (23a) all three methods
signal.pvs.list = onething(1)
