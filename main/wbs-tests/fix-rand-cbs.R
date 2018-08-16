## Synopsis: See randomized tests power compared to the nonrandomized version
library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")

## dosim_compare_old = dosim_compare ## this uses the double number of steps for CBS

nsim=500
oldresults = mclapply(1:nsim,function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    visc = c(5)
    results = dosim_compare_old(type="cbs.rand", n=30, lev=0, numIS=40, visc=visc,numSteps=4)
    return(results)
}, mc.cores=4)

nsim=500
newresults = mclapply(1:nsim,function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    visc = c(5)
    results = dosim_compare(type="cbs.rand", n=10, lev=0, numIS=40, visc=visc, numSteps=4)
    return(results)
}, mc.cores=4)



## allresults = lapply(allresults, function(myresult){names(myresult)=c("pvs","locs"); return(myresult)})
oldresults = lapply(oldresults, function(myresult){names(myresult)=c("pvs","locs"); return(myresult)})
newresults = lapply(newresults, function(myresult){names(myresult)=c("pvs","locs"); return(myresult)})
tab = do.call(rbind,allresults)
tab = tab[which(apply(tab,1,function(myrow)!all(is.na(myrow)))),]
## visc = n/2 + c((-5):5)
visc = n/2 + c((-1):1)
qqunif(tab[which(abs(tab[,"loc.wbs"]) %in% visc), "pvs"])
qqunif(tab[, "pvs"])

a=(tab[which(abs(tab[,"loc.wbs"]) %in% visc), "pvs"])


