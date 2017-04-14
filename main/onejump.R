## One-jump case.
onesim <- function(isim, sigma){

    ## generate data
    y <- c(rep(0,10),rep(1,10)) + rnorm(20,0,sigma)
    n = 20
    
    ###########################
    ## Do SBS-FS inference ####
    ###########################
    numSteps = 1
    method <- binSeg_fixedSteps
    obj <- method(y, numSteps)
    poly <- polyhedra(obj)
    p.bsfs = rep(NA,length(obj$cp))
    names(p.bsfs) = obj$cp
    contrast = list()
    for(ii in 1:length(obj$cp)){
        contrast[[ii]] <- contrast_vector(obj, ii)
        p.bsfs[ii] <- poly.pval(y=y, G=poly$gamma, u=poly$u, v=contrast[[ii]],sigma=sigma, bits=100)$pv
    }
    
    
    ###########################
    ## Do WBS-FS inference ####
    ###########################
    numSteps=1
    nsim.is = 100
    method <- wildBinSeg_fixedSteps
    numIntervals = 10
    intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)
    obj <- method(y, numSteps=numSteps, intervals=intervals) 
    poly <- polyhedra(obj)
    p.wbsfs = rep(NA,length(obj$cp))
    contrast <- make_all_segment_contrasts(obj)
    for(ii in 1:length(obj$cp)){
        p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
                                                v=contrast[[ii]], sigma=sigma,
                                                numSteps=numSteps,
                                                numIntervals=numIntervals,
                                                nsim.is=nsim.is, bits=100)
    }
    p.bsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.bsfs)
    p.wbsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.wbsfs)
    colnames(p.bsfs) = colnames(p.wbsfs) = c("isim","cp","pv")
    return(list(p.bsfs=p.bsfs, p.wbsfs=p.wbsfs))
}


## For a single sigma,
## 1. see if wbs randomized is more powerful.
n.std=5
std=seq(from=0,to=3,length=n.std)
nsim=100
results <- list(n.std)
for(i.std  in 1:n.std){
    manysimresult = lapply(1:nsim,function(isim){print(isim);onesim(isim,std[i.std])})
    plist.bsfs <- lapply(manysimresult, function(a)a$p.bsfs)
    plist.wbsfs <- lapply(manysimresult, function(a)a$p.wbsfs)
    pmat.bsfs = do.call(rbind, plist.bsfs)
    pmat.wbsfs = do.call(rbind, plist.wbsfs)
    results[[i.std]] <- list(pmat.bsfs, pmat.wbsfs)
}

save(results, file = "../data/onejump.Rdata")
