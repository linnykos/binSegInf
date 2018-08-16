### Synopsis: See recovery property and inference comparison of wbs vs cbs vs bs (hence (X)BS)
library(genlassoinf)
outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")
source("../main/wbs-tests/xbs-recovery-helpers.R")

### Run the actual simulation
levs=seq(from=0,to=4, by=.5)
nsim=1000
results.list = lapply(levs, function(mylev){
    print(mylev)
    a = mclapply(1:nsim,function(isim){
        printprogress(isim,nsim);
        onesim_recovery(mylev)},
        mc.cores=3)
    return(a)
})

filename="recovery-wbs-vs-cbs-vs-bs.Rdata"
save(list=c("results.list","levs"),
     file=file.path(outputdir,filename))

### Load and plot
load(file=file.path(outputdir,filename))
clean <- function(a){
    tab = do.call(rbind,a)
    apply(tab, 2, mean, na.rm=TRUE)
}
errortab = sapply(results.list, clean)
pdf(file.path(outputdir, "recovery-wbs-vs-cbs-vs-bs.pdf"), width=5, height=5)
lty = c(1,1,1)
cols = RColorBrewer::brewer.pal(3,"Set2")[c(1,2,3)]
matplot(x=levs, y=t(errortab[c(1,3,5),]), type='l', col=cols, lwd=2, lty=lty,
        ylim=c(0,1), axes=FALSE, ylab="", xlab=expression(delta))
axis(1)
axis(2)

legend("right", lty=lty, col=cols, legend=c("BS detection %",
                                            "WBS detection %",
                                            "CBS detection %"))
graphics.off()

