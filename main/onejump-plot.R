load(file="../data/onejump.Rdata")

## Plot the p-values 
myextract <- function(myresult, objname, locs){
    mysubresult = myresult[[objname]]
    inds = which((mysubresult)[,"cp"]%in%locs)
    return(mysubresult[inds,"pv"])
}


par(mfrow=c(2,3))
for(ii in 1:length(results)){
    qqunif(myextract(myresult = results[[ii]],
                     objname = "pmat.bsfs",
                     locs = n/2), plot.it = TRUE, ylim=c(0,1),xlim=c(0,1))
    xy = qqunif(myextract(myresult = results[[ii]],
                     objname = "pmat.wbsfs",
                     locs = n/2), xlim =c(0,1), plot.it=FALSE)
    points(xy$y~xy$x, col='red')
    title(main=paste("jump size =", levs[ii], "\n location=n/2"),ylim=c(0,1))
    legend("bottomright", col=c("black","red"), pch = c(1,1), legend = c("Binseg","WildBinSeg"))
}


