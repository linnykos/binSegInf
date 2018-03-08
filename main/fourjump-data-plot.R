## Synopsis: Plotting a four jump data example
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")
n=200
set.seed(0)
sigma = 1
y0    = fourjump(n=n, lev=2) + rnorm(n,0,1)
beta0 = lapply(c(0:2), function(lev) fourjump(lev=lev, n=n))

xlab = "Location"
w = 6; h = 4
ylim = c(-6,4)
mar = c(3.5,3.5,0.5,0.5)
xlim = c(0,n+10)

xticks = seq(from=0,to=n, by=50)
ltys.sig = c(2,2,1)
lwd.sig = 2
pch.dat = 16
pcol.dat = "grey50"
pch.contrast = 17
lty.contrast = 2
lcol.sig = 'red'
pcol.spike=3
pcol.segment=4
pcols.delta =   pcols.oneoff = RColorBrewer::brewer.pal(n=3,name="Set2")
pch.spike = 15
pch.segment = 17
cex.contrast = 1.2
cex.dat = .8

pdf(file.path(outputdir,"fourjump-data.pdf"), width=w, height=h)
par(mar=c(4.1,3.1,3.1,1.1))
plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch.dat, col=pcol.dat, cex=cex.dat);
axis(1, at = xticks, labels = xticks); axis(2)
for(ii in 1:3) lines(beta0[[ii]],col="red",lty=ltys.sig[ii], lwd=lwd.sig)
for(ii in 0:2) text(x=30,y=ii+.2, label = bquote(delta==.(ii)))
legend("bottomleft", pch=c(pch.dat,NA),
       lty=c(NA,1), lwd=c(NA,2),
       col = c(pcol.dat, lcol.sig),
       pt.cex = c(cex.contrast, NA),
       legend=c("Data", "Mean"))
## title(main=expression("Data example"))
graphics.off()
