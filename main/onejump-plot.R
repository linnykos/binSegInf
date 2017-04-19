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
    qqunif(myextract(myresult = results[[ii]],
                     objname = "pmat.bsfs",
                     locs = sim.settings$n/2), plot.it = TRUE, ylim=c(0,1),xlim=c(0,1))
    xy = qqunif(myextract(myresult = results[[ii]],
                     objname = "pmat.wbsfs",
                     locs = sim.settings$n/2), xlim =c(0,1), plot.it=FALSE)
    points(xy$y~xy$x, col='red')
    title(main=paste("jump size =", sim.settings$levs[ii], "\n location=n/2"),ylim=c(0,1))
    legend("bottomright", col=c("black","red"), pch = c(1,1), legend = c("Binseg","WildBinSeg"))
}



#########################
### Artificial signal ###
#########################

## load(file="~/Desktop/data/artificial.Rdata")
load(file="../data/artificial-smaller.Rdata")

newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
mndf = makeDmat(length(newmn2),type="tf",ord=0)%*%newmn2
mndf[abs(mndf)<1E-10]=0
cp <- which(mndf!=0)




par(mfrow=c(4,4))
n = sim.settings$n
for(ii in 1:length(results)){

  for(jj in 1:length(cp)){
    qqunif(myextract(myresult = results[[ii]],
                     objname = "pmat.bsfs",
                     locs = cp[jj]), plot.it = TRUE, ylim=c(0,1),xlim=c(0,1))
  
    xy = qqunif(myextract(myresult = results[[ii]],
                     objname = "pmat.wbsfs",
                     locs = cp[jj]), xlim =c(0,1), plot.it=FALSE)
    points(xy$y~xy$x, col='red')
    title(main=paste("jump size =", sim.settings$levs[ii], "\n location=", cp[jj]),ylim=c(0,1))
    legend("bottomright", col=c("black","red"), pch = c(1,1), legend = c("Binseg","WildBinSeg"))
}
}
