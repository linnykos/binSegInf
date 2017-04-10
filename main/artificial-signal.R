library(DNACopy)
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
mn2 = softthresh(a,lambda=cv$lambda.1se,gamma=0.1)

## Plot
sigma=.01
set.seed(0)
y = mn2 + rnorm(length(mn2),0, sigma)
plot(y)

## Do SBS-FS inference
method <- binSeg_fixedSteps
obj <- method(y, 6)
poly <- polyhedra(obj)
p = rep(NA,length(obj$cp))
names(p) = obj$cp
for(ii in 1:length(obj$cp)){
    contrast <- contrast_vector(obj, ii)
    p[ii] <- pvalue(y, poly, contrast)
}
plot(y)
abline(v=names(p))
text(x=names(p),y=rep(0.4 + rnorm(length(p),0,0.1),length(p)), label=(round(p,3)))
##  D





## method <- wildBinSeg_fixedThresh
## obj <- method(y, thresh=0.15, numIntervals=100)
