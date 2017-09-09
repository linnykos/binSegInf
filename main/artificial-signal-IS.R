load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    ## newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
## n = length(coriell_mn(1))
n = 100
numSteps = 1
numIntervals = 100
lev=3
sigma=1
nreplicate =  c(100,50,20)
nsims.is = c(10,100,500)
mc.cores=5
sim.settings = list(numIntervals=numIntervals,
                    nreplicate=nreplicate,
                    lev=lev,
                    numSteps=numSteps,
                    nsim.is=nsim.is)

## ## Temporary settings
## nreplicate =  c(100,50,10)/10
## nsims.is=c(3,3,3)
## sigma=1
## n=20


for(ii in 1:3){
    nsim.is = nsims.is[ii]
    cat("nsim.is is", nsim.is, "out of", nsims.is, fill=TRUE)

    ## Generate data
    set.seed(0)
    ## y <- coriell_mn(lev,n) + rnorm(n,0,std)
    y <- mn.onejump(lev,n) + rnorm(n,0,1)

    ## p.wbsfs.list = list()
    ## for(irep in 1:(nreplicate[ii])){
    p.wbsfs.list = mclapply(1:(nreplicate[ii]), function(irep){

        cat("\r", "replicate", irep, "out of", nreplicate[ii])

        method <- wildBinSeg_fixedSteps
        intervals <- generate_intervals(n=length(y),numIntervals=numIntervals, seed=0)
        obj <- method(y, numSteps=numSteps, intervals=intervals)
        poly <- polyhedra(obj)

        contrast <- make_all_segment_contrasts(obj)
        p.wbsfs = rep(NA,length(obj$cp))
        names(p.wbsfs) = obj$cp
        for(ii in 1:length(obj$cp)){
            p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
                                                    v=contrast[[ii]], sigma=sigma,
                                                    numSteps=numSteps,
                                                    numIntervals=numIntervals,
                                                    nsim.is=nsim.is, bits=100)
        }
        ## p.wbsfs.list[[irep]] = p.wbsfs
    }, mc.cores=mc.cores)

    ## Save each result
    filename = paste0("artificial-nsimis-", nsim.is ,".Rdata")
    save(p.wbsfs.list, sim.settings, file = file.path("../results", filename) )
    cat(fill=TRUE)
}


## ## In a single run, the polyhedron is too big. memory won't allow it to run.
## ## Find out why (testing code should be standalone.)
## coriell_mn <- function(lev=1,n){
##     h = max(abs(newmn))
##     return((newmn / h * std) * lev)
## }

## ptm <- proc.time()
## n = length(coriell_mn(1))
## set.seed(0)
## lev=3
## numIntervals=1000
## numSteps=1
## augment=TRUE
## n = length(newmn)
## y <- coriell_mn(lev,n) + rnorm(n,0,std)
## p.wsbs.list = list()
## method <- wildBinSeg_fixedSteps
## intervals <- generate_intervals(n=length(y),numIntervals=numIntervals, seed=0)
## obj <- method(y, numSteps= numSteps, intervals=intervals, augment=augment)

## print("ran method")
## print(proc.time() - ptm)

## poly <- polyhedra(obj)

## print("collected polyhedron")
## print(proc.time() - ptm)

## contrast <- make_all_segment_contrasts(obj)
## ## for(ii in 1:length(obj$cp)){
## ii==1
## p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
##                                         v=contrast[[ii]], sigma=sigma,
##                                         numSteps=numSteps,
##                                         numIntervals=numIntervals,
##                                         nsim.is=nsim.is, bits=100)
