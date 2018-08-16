## Synopsis: Run a single example on original 05296 data.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir, filename))
source(file=file.path("../main/artificial/artif-helpers.R"))


## Remaining data after using first 200 for noise estimation
sigma = sd(y.orig[1:200])
y = y.orig[201:length(y.orig)]

## Plot settings
load(file.path(outputdir,"rwbs-orig-segment.Rdata"))
w = 10; h = 4
xlab = "Location"
ylab = ""
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-1, 1.3)
mar = c(4.5,4.5,0.5,0.5)
obj = orig.pvs.rwbs.segment
cp = obj$locs.all
Letters = toupper(letters[1:length(cp)])
segments = lapply(1:(length(cp)+1), function(ii){
    v = c(0,sort(abs(cp)),length(y));
    (v[ii]+1):(v[ii+1]) })
segment.means = sapply(segments, function(mysegment){mean(y[mysegment])})
mn0 = do.call(c,sapply(1:length(segments), function(ii){rep(segment.means[ii], length(segments[[ii]]))}))


## Conduct naive tests and collect them.
pval.naive.ryan = function(y0, loc, loc.left, loc.right, sigma) {
  gap = abs(mean(y0[loc.left:(loc-1)]) - mean(y0[loc:(loc.right-1)]))
  sd = sqrt(sigma^2*(1/(loc-loc.left)+1/(loc.right-loc)))
  return(2*(1-pnorm(gap,mean=0,sd=sd)))
}
abs.cp= sort(abs(cp))
aug.abs.cp = c(0, abs.cp, length(y))
lefts = aug.abs.cp[-c(length(aug.abs.cp), length(aug.abs.cp)-1)]
mids = abs.cp
rights = aug.abs.cp[-c(1:2)]
naive.pvs = Map(function(left,mid,right){pval.naive.ryan(y,mid,left,right,sigma)}, lefts, mids, rights)
names(naive.pvs) = cp[order(abs(cp))]
unlist(naive.pvs)



## Make the plot
pdf(file.path(outputdir,"05296-intro.pdf"), width=w, height=h)
par(mar=mar)
plot(y,ylim=ylim,xlab=xlab, col=pcol,pch=pch,axes=FALSE, ylab=ylab)
axis(1); axis(2)
lines(mn0, col='red', lwd=2)
abline(v=abs(cp),col='lightgrey')
rect(xleft=0, xright=length(y.orig), ybottom=1.2, ytop=2, border=NA, col="white")
text(x=abs(cp)[-c(9,4)], y=1.2*max(y),label=Letters[-c(9,4)],cex=1)
text(x=abs(cp)[9]+10, y=1.2*max(y),label=Letters[9],cex=1)
text(x=abs(cp)[4]+10, y=1.2*max(y),label=Letters[4],cex=1)
legend("bottomright", pch=c(pch,NA,NA), lty=c(NA,1,2),lwd=c(NA,lwd,1),
       col = c(pcol,"red","black"),legend=c("aCGH Data","Estimate","Changepoint"), bg="white")
graphics.off()

## Make the table
obj = all.results[[filenames[1]]]
pvs = obj$pvs
pvs = rbind(round(pvs,3))
colnames(pvs) = Letters[order(abs(obj$locs.all))]
pvs = rbind(naive.pvs, pvs)
rownames(pvs) = c("naive", "adjusted")
xtable::xtable(pvs)


