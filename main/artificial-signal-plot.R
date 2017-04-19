load("../data/artificial.Rdata")

## Extract results
results[[1]]$pmat.bsfs[,"cp"]
results[[1]]$pmat.wbsfs[,"cp"]
results[[1]]$pmat.wbsfs.nonrand[,"cp"]

## Example plot
cp = results[[2]]$pmat.bsfs[,"cp"]
pv = results[[2]]$pmat.bsfs[,"pv"]
mn0 = sim.settings$mn(lev=sim.settings$levs[2],n=sim.settings$n)
y0 = mn0 + rnorm(sim.settings$n, 0, sim.settings$sigma)
plot(y0)
lines(mn0, col='red')
abline(v=cp,col='lightgrey', 3)
text(x=cp, y = rep(0.3, length(cp))+runif(length(cp),0,0.1),label=round(pv,3),cex=2)
