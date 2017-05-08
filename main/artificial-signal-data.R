library(DNAcopy)
library(genlasso)
library(genlassoinf)
data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23
CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))
y.orig = (coriell[,4])
y.orig = y.orig[!is.na(y.orig)]
a = fusedlasso1d(y.orig)
cv = cv.trendfilter(a)

## Harvest a sparse-fused lasso
mn.elnet = softthresh(a,lambda=cv$lambda.1se,gamma=0.05)

## Refit the segment means (because sparse fused lasso estimate is biased)
D = makeDmat(length(y.orig),type="tf",ord=0)
ydiff = D %*% mn.elnet
ydiff[abs(ydiff)<1E-10]=0
cp <- which(ydiff!=0)
cp <- cp[c(1,4,6,7,8)]
## cleanmn <- get_piecewise_mean(y,cp)
## std = sd(y-cleanmn)
## segment.means = sapply(segments, function(mysegment){mean(y[mysegment])})
segments = lapply(1:(length(cp)+1), function(ii){ v=c(0,cp,length(y.orig));(v[ii]+1):(v[ii+1]) })
segment.means = sapply(segments, function(mysegment){mean(y.orig[mysegment])})

## Make segment means
cleanmn = rep(NA,length(y.orig))
lapply(1:length(segments), function(ii){cleanmn[segments[[ii]]] <<- segment.means[ii]})
std = sd(y.orig-cleanmn)


## Make means, forcing some things to be zero.
newmn = rep(NA,length(y.orig))
segment.means[c(1,3,5)] = 0
lapply(1:length(segments), function(ii){newmn[segments[[ii]]] <<- segment.means[ii]})
newmn.short = newmn[1101:1300]
newmn.shorter = newmn.short[seq(from=1,to=200,length=100)]

## Save it
resid.cleanmn <- y.orig - cleanmn
filename = "../data/coriell.Rdata"
save(y.orig,
     segment.means,
     segments,
     cleanmn,
     resid.cleanmn,
     newmn,
     newmn.short,
     std,
     file = filename)
