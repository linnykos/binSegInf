load("../data/artificial-smaller.Rdata")

## Extract results
results[[1]]$pmat.bsfs[,"cp"]
results[[1]]$pmat.wbsfs[,"cp"]
results[[1]]$pmat.wbsfs.nonrand[,"cp"]


## Plot settings
xlab = "Location"
ylab = ""
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-1.5,2)
mar = c(4.5,4.5,0.5,0.5)
cp = sort(results[[ii]]$pmat.bsfs[,"cp"])
beta0 <- get_piecewise_mean(y, cp)
Letters = toupper(letters[1:length(cp)])

## Example plot
ylim = c(range(y0*1.5))
cp = results[[ii]]$pmat.bsfs[,"cp"]
pv = results[[ii]]$pmat.bsfs[,"pv"]
mn0 = sim.settings$mn(lev=sim.settings$levs[ii],n=sim.settings$n)
y0 = mn0 + rnorm(sim.settings$n, 0, sim.settings$sigma)
## pdf("artif1.pdf", width=w,height=h)
par(mar=mar)
plot(y0,ylim=ylim,xlab=xlab, col=pcol,pch=pch,axes=FALSE, ylab=ylab)
axis(1);axis(2)
lines(mn0, col='red')
abline(v=cp,col='lightgrey')
text(x=cp, y=1.2*max(y0),label=Letters,cex=1)
round(pv,3)
rep(max(y0)*.8, length(cp))+runif(length(cp),0,0.1)
legend("bottomright", pch=c(pch,NA,NA), lty=c(NA,1,2),lwd=c(NA,lwd,1),
       col = c(pcol,"blue","black"),legend=c("aCGH Data","Estimate","Changepoint"), bg="white")
## graphics.off()





## #########################
## ### Artificial signal ###
## #########################

## ## load(file="~/Desktop/data/artificial.Rdata")
## load(file="../data/artificial-smaller.Rdata")

## newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
## mndf = makeDmat(length(newmn2),type="tf",ord=0)%*%newmn2
## mndf[abs(mndf)<1E-10]=0
## cp <- which(mndf!=0)




## par(mfrow=c(4,4))
## n = sim.settings$n
## for(ii in 1:length(results)){

##   for(jj in 1:length(cp)){
##     qqunif(myextract(myresult = results[[ii]],
##                      objname = "pmat.bsfs",
##                      locs = cp[jj]), plot.it = TRUE, ylim=c(0,1),xlim=c(0,1))

##     xy = qqunif(myextract(myresult = results[[ii]],
##                      objname = "pmat.wbsfs",
##                      locs = cp[jj]), xlim =c(0,1), plot.it=FALSE)
##     points(xy$y~xy$x, col='red')
##     title(main=paste("jump size =", sim.settings$levs[ii], "\n location=", cp[jj]),ylim=c(0,1))
##     legend("bottomright", col=c("black","red"), pch = c(1,1), legend = c("Binseg","WildBinSeg"))
## }
## }
