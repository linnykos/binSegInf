load("../data/artificial.Rdata")

## Extract results
results[[1]]$pmat.bsfs[,"cp"]
results[[1]]$pmat.wbsfs[,"cp"]
results[[1]]$pmat.wbsfs.nonrand[,"cp"]

## Example plot
cp = results[[1]]$pmat.bsfs[,"cp"]
pv = results[[1]]$pmat.bsfs[,"pv"]
plot(y)
lines(newmn)
abline(v=cp,col='red')
text(x=cp, y = rep(0.3, length(cp))+runif(length(cp),0,0.1),label=round(pv,3),cex=2)
