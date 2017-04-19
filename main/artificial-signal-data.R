library(DNAcopy)
library(genlasso)
library(genlassoinf)
data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23
CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

y = (coriell[,4])
y = y[!is.na(y)]
a = fusedlasso1d(y)
cv = cv.trendfilter(a)

## Harvest a sparse-fused lasso
mn.elnet = softthresh(a,lambda=cv$lambda.1se,gamma=0.05)

## Refit the segment means (because sparse fused lasso estimate is biased)
D = makeDmat(length(y),type="tf",ord=0)
ydiff = D %*% mn.elnet
ydiff[abs(ydiff)<1E-10]=0
cp <- which(ydiff!=0)
cp <- cp[c(1,4,6,7,8)]
segments = lapply(1:(length(cp)+1), function(ii){ v=c(0,cp,length(y));(v[ii]+1):(v[ii+1]) })
segment.means = sapply(segments, function(mysegment){mean(y[mysegment])})

## Make segment means
cleanmn = rep(NA,length(y))
lapply(1:length(segments), function(ii){cleanmn[segments[[ii]]] <<- segment.means[ii]})
std = sd(y-cleanmn)

## Make means, forcing some things to be zero.
newmn = rep(NA,length(y))
segment.means[c(1,3,5)] = 0
lapply(1:length(segments), function(ii){newmn[segments[[ii]]] <<- segment.means[ii]})
newmn.short = newmn[1001:1500]

## Save it
resid.cleanmn <- y - cleanmn
filename = "../data/coriell.Rdata"
y.orig<-y
save(y.orig,
     segment.means,
     segments,
     cleanmn,
     resid.cleanmn,
     newmn,
     newmn.short,
     std,
     file = filename)
