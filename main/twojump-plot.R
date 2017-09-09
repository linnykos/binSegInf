load(file="../results/twojump.Rdata")

## Plot the p-values
myextract <- function(myresult, objname, locs){
    mysubresult = myresult[[objname]]
    inds = which((mysubresult)[,"cp"]%in%locs)
    return(mysubresult[inds,"pv"])
}


par(mfrow=c(2,3))
n = sim.settings$n
## Plot settings
mar = c(4.5,4.5,0.5,0.5)
for(ii in 1:length(results)){
    pmat.bsfs = results[[ii]]$pmat.bsfs
    pmat.wbsfs = results[[ii]]$pmat.wbsfs
    pmat.wbsfs.plain = results[[ii]]$pmat.wbsfs.plain
    loc=2*(sim.settings$n)/3
    pdf(file.path("../main/figures", paste0("twojump-",ii,".pdf")),width=5,height=5)
    par(mar=mar)
    qqunif(pmat.bsfs[,loc],ylim=c(0,1),xlim=c(0,1))
    qqunif_add(pmat.wbsfs[,loc],col='red',ylim=c(0,1),xlim=c(0,1))
    qqunif_add(pmat.wbsfs.plain[,loc],col='green',ylim=c(0,1),xlim=c(0,1))
    title(main=paste("jump size =", sim.settings$levs[ii], "\n location=",loc),
          ylim=c(0,1))
    legend("bottomright", col=c("black","red", "green"), pch = c(1,1,1),
           legend = c("Binseg","WildBinSeg", "WildBinSeg-plain"))
    graphics.off()
}


