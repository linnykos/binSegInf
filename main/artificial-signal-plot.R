load("../data/artificial-smaller.Rdata")

## Extract results
results[[1]]$pmat.bsfs[,"cp"]
results[[1]]$pmat.wbsfs[,"cp"]
results[[1]]$pmat.wbsfs.nonrand[,"cp"]


ii=2

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
