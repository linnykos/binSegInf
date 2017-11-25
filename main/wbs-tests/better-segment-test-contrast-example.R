## Synopsis: Make a single example showing paper and contrasts for the paper

source("../main/wbs-tests/sim-helpers.R")
n=200
beta0 = fourjump(lev=1,n=n)

xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
ylim = c(-5,7)
mar = c(4.5,4.5,0.5,0.5)
xlim = c(0,n+10)
let = c("A","B")
ltys.sig = c(2,2,1)
lwd.sig = 2
pch.dat = 16
pcol.dat = "grey80"
pch.contrast = 17
lty.contrast = 2
lcol.sig = 'red'
pcol.orig=3
pcol.plus=4
pcols.delta =   pcols.oneoff = RColorBrewer::brewer.pal(n=3,name="Set2")
pch.orig = 15
pch.plus = 17
cex.contrast = 1.2

pdf(file.path(outputdir,"better-segment-example.pdf"), width=5,height=5)

## Run Wild binary segmentation for four steps
g = wildBinSeg_fixedSteps(y=y0, numSteps=4, numIntervals=n, inference.type="none")
locs = c(120+c((-5):5))
## locs = c(80+c((-5):5))
ii = which(g$cp %in% locs)
if(length(ii)!=1)  next
test.cp = g$cp[which(g$cp %in% locs)]
test.cp.sign = g$cp.sign[which(g$cp %in% locs)]
wbs.ends = g$results[ii,c("max.s", "max.e")]

## Get contrasts
which.closest = order(abs(g$cp - g$cp[ii]))[2:3]
closest.ends = sort(c(0,g$cp,n)[order(abs(c(0,g$cp,n) - g$cp[ii]))[2:3]])
vlist = make_all_segment_contrasts_from_cp(cp=g$cp,n=n, cp.sign=g$cp.sign)
vlist2 = make_all_segment_contrasts_from_wbs(wbs_obj=g)
v = vlist[[toString(test.cp*test.cp.sign)]]
v.wbs = vlist2[[toString(test.cp*test.cp.sign)]]
v.wbs[which(v.wbs==0)]=NA
v[which(v==0)]=NA
fac=30
v = v*fac
v.wbs = v.wbs*fac
xcoord=1:n

## Plot them
par(mar=c(4.1,3.1,3.1,1.1))
plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol.dat);
axis(1);axis(2)
lines(beta0,col="red",lty=ltys.sig[3], lwd=lwd.sig)
points(v~xcoord, pch = pch.orig, col = pcol.orig)
points(v.wbs~xcoord, pch = pch.plus, col = pcol.plus)
abline(v=g$cp[ii],col='lightgrey')
legend("topleft", pch=c(pch.dat,NA,pch.orig,pch.plus), lty=c(NA,1,NA,NA),
       lwd=c(NA,2,NA,NA), col = c(pcol.dat, lcol.sig, pcol.orig, pcol.plus),
       pt.cex = c(cex.contrast, NA, cex.contrast, cex.contrast),
       legend=c("Data", "Mean","Segment contrast", "Segment+ contrast"))
## title(main=expression("Data example"))


## Add arrows for signal size
offset=4
y.offset=0.2
arrows(x0=40-offset, y0=0, x1=40-offset, y1 = 1, col="black", length=.1, angle=20, lwd=2,code=2)
text(x=40-offset*2, y=1+y.offset, label=expression(delta==1))
arrows(x0=160+offset, y0=0, x1=160+offset, y1 = -2, col="black", length=.1, angle=20, lwd=2,code=2)
text(x=160+offset*2, y=-2-y.offset, label=expression(2~delta==2))


graphics.off()
