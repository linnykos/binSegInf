## Synopsis: prepare artificial mean for Coriell Cell lines 05296 and 13330.
## Both are used in the original CBS paper, available from DNACopy R package,
## and have known truths.
library(DNAcopy)
library(genlasso); library(genlassoinf)
data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23
CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))
names(coriell)
y.orig = (coriell[,"Coriell.05296"])

## Harvest a sparse-fused lasso
nonmissing = which(!is.na(y.orig))
y.orig = y.orig[nonmissing]
coriell = coriell[nonmissing,]
a = fusedlasso1d(y.orig)
cv = cv.trendfilter(a)
mn.elnet = softthresh(a,lambda=cv$lambda.1se,gamma=0.05)
plot(mn.elnet)

## Refit the segment means (because sparse fused lasso estimate is biased)
D = makeDmat(length(y.orig),type="tf",ord=0)
ydiff = D %*% mn.elnet
ydiff[abs(ydiff)<1E-10]=0
cp <- which(ydiff!=0)
cp <- cp[c(1,4,6,7,8)]
segments = lapply(1:(length(cp)+1), function(ii){
    v = c(0,cp,length(y.orig));
    (v[ii]+1):(v[ii+1]) })
segment.means = sapply(segments, function(mysegment){mean(y.orig[mysegment])})
cleanmn = rep(NA,length(y.orig))
for(ii in 1:length(segments){
    cleanmn[segments[[ii]]] <<- segment.means[ii]
}
std = sd(y.orig-cleanmn)

## Make means, forcing some things to be zero.
newmn = rep(NA,length(y.orig))
segment.means[c(1,3,5)] = 0
lapply(1:length(segments), function(ii){newmn[segments[[ii]]] <<- segment.means[ii]})
resid.cleanmn <- y.orig - cleanmn

## Save it
filename = "../data/coriell05296.Rdata"
## filename = "../data/coriell13330.Rdata"
save(y.orig,
     segment.means,
     segments,
     cleanmn,
     resid.cleanmn,
     newmn,
     std,
     file = filename)



## plot(y.orig)
## y.orig = (coriell[,"Coriell.13330"])
## plot(y.orig)
## dat = y.orig[1000:2000]
## dat = (dat-mean(dat))/sd(dat)
## qqnorm(dat)
## abline(0,1)




## ## Plot the data once.
## plot(y.orig)
## boundary.inds = c(0,cumsum(sapply(1:23, function(ii)sum(coriell[,"Chromosome"]==ii))))
## xt = sapply(1:(length(boundary.inds)-1), function(ii)mean(boundary.inds[c(ii,ii+1)] ))
## abline(v=boundary.inds)
## text(x=xt,y=rep(1,23), label=1:23)
