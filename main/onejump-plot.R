load(file="~/Desktop/data/onejump.Rdata")

## Plot the p-values
myextract <- function(myresult, objname, locs){
    mysubresult = myresult[[objname]]
    inds = which((mysubresult)[,"cp"]%in%locs)
    return(mysubresult[inds,"pv"])
}


par(mfrow=c(2,3))
n = sim.settings$n
for(ii in 1:length(results)){
    pmat.bsfs = results[[1]]$pmat.bsfs
    pmat.wbsfs = results[[1]]$pmat.wbsfs
    pmat.wbsfs.plain = results[[1]]$pmat.wbsfs.plain
    loc=(sim.settings$n)/2
    qqunif(pmat.bsfs[,loc])
    qqunif_add(pmat.wbsfs[,loc],col='red')
    qqunif_add(pmat.wbsfs.plain[,loc],col='green')
    title(main=paste("jump size =", sim.settings$levs[ii], "\n location=n/2"),
          ylim=c(0,1))
    legend("bottomright", col=c("black","red", "green"), pch = c(1,1,1),
           legend = c("Binseg","WildBinSeg", "WildBinSeg-plain"))
}


