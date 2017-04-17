library(DNAcopy)
library(genlasso)
data(coriell)
#Combine into one CNA object to prepare for analysis on Chromosomes 1-23
CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

head(coriell)
plot(coriell[coriell[1,]==1,4])
plot(coriell[,4])

y = (coriell[,4])
y = y[!is.na(y)]
a = fusedlasso1d(y)
plot(a)
plot(y)

## Choose lambda by 5-fold cross-validation
cv = cv.trendfilter(a)
plot(cv)
plot(a,lambda=cv$lambda.1se,main="One standard error rule")
mn = as.numeric((coef(a,lambda=cv$lambda.1se))$beta)

## Harvest a sparse-fused lasso
## mn2 = softthresh(a,lambda=cv$lambda.1se,gamma=0.1)
mn2 = softthresh(a,lambda=cv$lambda.1se,gamma=0.05)
## mn2 = c(rep(0,100), rep(1,100))

## Refit the segment means (because sparse fused lasso estimate is biased)
D = makeDmat(length(y),type="tf",ord=0)
ydiff = D %*% mn2
ydiff[abs(ydiff)<1E-10]=0

## The mean doesn't look very flattering anymore. But we proceed
plot(y)
lines(mn2,col='red')
abline(v=which(ydiff!=0))
abline(v=seq(from=0,to=2000,by=100),col='lightgrey')
cp <- which(ydiff!=0)
cp <- cp[c(1,3,4,6,7,8)]
segments = lapply(1:(length(cp)+1), function(ii){ v=c(0,cp,length(y));(v[ii]+1):(v[ii+1]) })
newmn = rep(NA,length(y))
lapply(segments, function(mysegment){newmn[mysegment] <<- mean(y[mysegment])})


## Plot
sigma =sd(y-mn)
set.seed(0)
y = newmn + rnorm(length(newmn),0, sigma)
plot(y)

## Save it
mn <- newmn
resid <- y-mn
filename = "data/coriell.Rdata"
y.orig=y
save(y.orig, mn, resid, file = filename)
