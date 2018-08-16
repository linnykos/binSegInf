## Synopsis: Segment+ test endpoints are simulated and scatterplots are produced

## Run simulation to gather endpoints
source("../main/wbs-tests/sim-helpers.R")
sigma = 1
n=200
nsim=500
pvmat = matrix(ncol=2,nrow=nsim)
colnames(pvmat) = c("orig", "wbs")
endsmat = matrix(ncol=5, nrow=nsim)
colnames(endsmat) = c("bk", "wbs.l", "wbs.r", "ord.l", "ord.r")
for(isim in 274:nsim){
    printprogress(isim,nsim)
    set.seed(isim)
    y0 = fourjump(lev=1,n=n) + rnorm(n,0,sigma)
    g = wildBinSeg_fixedSteps(y=y0, numSteps=4, numIntervals=n, inference.type="rows")

    locs = c(120+c((-5):5))
    ## locs = c(80+c((-5):5))
    ii = which(g$cp %in% locs)
    if(length(ii)!=1)  next
    test.cp = g$cp[which(g$cp %in% locs)]
    test.cp.sign = g$cp.sign[which(g$cp %in% locs)]
    wbs.ends = g$results[ii,c("max.s", "max.e")]

    which.closest = order(abs(g$cp - g$cp[ii]))[2:3]
    closest.ends = sort(c(0,g$cp,n)[order(abs(c(0,g$cp,n) - g$cp[ii]))[2:3]])

    ## Compare p-values
    vlist = make_all_segment_contrasts_from_cp(cp=g$cp,n=n, cp.sign=g$cp.sign)
    vlist2 = make_all_segment_contrasts_from_wbs(wbs_obj=g)
    v = vlist[[toString(test.cp*test.cp.sign)]]
    v.wbs = vlist2[[toString(test.cp*test.cp.sign)]]
    poly = polyhedra(obj=g$gamma,u=rep(0,nrow(g$gamma)))
    pv = poly.pval2(v=v,sigma=sigma,poly=poly,y=y0)$pv
    pv.wbs = poly.pval2(v=v.wbs,sigma=sigma,poly=poly,y=y0)$pv

    ## Store it
    pvmat[isim,] = c(pv,pv.wbs)
    endsmat[isim,] = c(g$cp[ii], wbs.ends, closest.ends)
}
filename = "better-segment-endpoints-120.Rdata"
## filename = "better-segment-endpoints-80.Rdata"
outputdir="../output"
save(list=c("pvmat", "endsmat"), file = file.path(outputdir, filename))
load(file = file.path(outputdir, filename))



## Plotting the ends via scatterplots
## filename = "better-segment-endpoints-80.Rdata"
filename = "better-segment-endpoints-120.Rdata"
outputdir="../output"
load(file = file.path(outputdir, filename))
w=5;h=5
## pdf(file=file.path(outputdir, "better-segment-endpoints-scatterplot-80.pdf"), width=w, height=h)
pdf(file=file.path(outputdir, "better-segment-endpoints-scatterplot-120.pdf"), width=w, height=h)
plot(endsmat[,"wbs.l"]~endsmat[,"ord.l"], ylab = "segment+", xlab="segment",
     axes=FALSE, ylim=c(0,n), xlim=c(0,n), pch=1)
points(endsmat[,"wbs.r"]~endsmat[,"ord.r"], col='red', pch=1)
legend('bottomright', col=c("black", "red"), pch=c(1,1), legend=c("left endpoints", "right endpoints"))
abline(v=c(40,80,120), col='lightgrey')
abline(h=c(40,80,120),col='lightgrey')
abline(0,1)
axis(1);axis(2)
graphics.off()


## See p-values
qqunif(pvmat[,"orig"], col='black', pch=16)
a = qqunif(pvmat[,"wbs"], plot.it=FALSE)
points(a, col='red', pch=16)
legend("bottomright", legend=c("orig", "improved"), col = c("black", "red"), pch=c(16,16))
