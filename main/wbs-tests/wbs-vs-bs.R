## Synopsis: See recovery property and inference comparison of wbs vs cbs vs bs
library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
source("../main/wbs-tests/wbs-vs-bs-helpers.R")

## Run the actual simulation
levs=seq(from=0,to=4, by=.5)
nsim=100
results.list = lapply(levs, function(mylev){
    print(mylev)
    a = mclapply(1:nsim,function(isim){printprogress(isim,nsim);onesim_recovery(mylev)}, mc.cores=4)
    return(a)
})
filename="recovery-wbs-vs-cbs-vs-bs.Rdata"
save(list=c("results.list","levs"),
     file=file.path(outputdir,filename))

## Load and plot
load(file=file.path(outputdir,filename))
clean<-function(a){
    tab = do.call(rbind,a)
    apply(tab, 2, mean, na.rm=TRUE)
}
errortab = sapply(results.list, clean)
pdf(file.path(outputdir, "recovery-wbs-vs-cbs-vs-bs.pdf"), width=5, height=5)
lty = c(1,2,1,2)
cols = RColorBrewer::brewer.pal(3,"Set2")[c(1,1,2,2)]
matplot(x=levs, y=t(errortab), type='l', col=cols, lwd=2, lty=lty, ylim=c(0,1))
legend("topleft", lty=lty, col=cols, legend=c("BS recovery", "BS error", "WBS recovery", "WBS error"))
graphics.off()


## Power comparison
levs=seq(from=0,to=4, by=.5)
nsim=100
results.list = lapply(levs, function(mylev){
    print(mylev)
    a = mclapply(1:nsim,function(isim){printprogress(isim,nsim);onesim_power(mylev)}, mc.cores=4)
    return(a)
})

source("../main/wbs-tests/wbs-vs-bs-helpers.R")
a = onesim_power(mylev)

