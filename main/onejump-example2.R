## Setup
outputdir = "output"
library(RColorBrewer)
library(igraph)
library(stringdist)
## library(BinSegInf)

## ## Load R package
rpackage_dir = "~/repos/binSegInf/binSegInf"
setwd(rpackage_dir)
library(devtools)
load_all()

## Source in the p-value functions from the selectiveInference package 
pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
}

## Setting 3
n = 20
sigma=1
lev1=0
lev2=3
thresh = 2
nsim=1000
pvals.list = list()
lambdas = seq(from=.5, to = 2, by = .25)
for(ilambda in 1:length(lambdas)){
    cat(ilambda, "out of", length(lambdas),"lambdas", '\n')
    lambda = lambdas[ilambda]
    pvals = matrix(NA, nrow=nsim, ncol=n)
    for(isim in 1:nsim){
        cat('\r', isim, "out of", nsim, "simulations")
        ## Set it up
        y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
        slist = elist = blist = Blist = zlist = Zlist =
            matrix(NA, nrow = n, ncol = 2^8)

        ## Do Binary Segmentation
        tryCatch({
            binseg(s = 1, e = n, j = 0, k = 1, thresh = lambda, y = y, n = n)}, error=function(e){})
    
        ## Bundle into one thing
        bs.output = list(slist = slist, elist = elist, blist = blist, Blist = Blist,
                         zlist = zlist, Zlist = Zlist, y = y, thresh = lambda)
    
        ## Get Gamma matrix and u vector
        Gu = get.polyhedron(bs.output, lambda, verbose=F)
    
        ## Conduct segment TG tests
        test.b.list = sort(collapse(blist))
        pvals[isim, test.b.list] = as.numeric(sapply(test.b.list,
                                          function(test.b){return( pval.fl1d(y = y,
                                                                             G = Gu$G,
                                                                             dik = make.v(test.b, bs.output),
                                                                             sigma = sigma,
                                                                             u = Gu$u))}))
    }
    pvals.list[[ilambda]] = pvals
}

sim.settings = list(lev1=lev1,lev2=lev2,thresh=thresh,nsim=nsim,sigma=sigma, n=n, lambdas=lambdas)
save(pvals.list, sim.settings, file='onejump-by-lambda.Rdata')



## Plotting code
load('onejump-by-lambda.Rdata')
pdf("~/Desktop/onejump-by-lambda.pdf")
par(mfrow=c(2,4))
lev2 = 3
set.seed(0)
y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
plot(y, pch = 16, main = "data and signal (example)")
lines(c(rep(lev1,n/2),rep(lev2,n/2)), col='red',lwd=2)
for(ilambda in 1:length(lambdas)){
hist(pvals.list[[ilambda]][,10],main="")
title(main=bquote(lambda==.(lambdas[ilambda])))
}
graphics.off()
## unif.sample = seq(from=0,to=1,length=nsim)
## qqplot(pvals.list[[1]][,10], unif.sample,col=mycol)

## replicate(100, do.once())

## If I want R to find stuff well, read this:
## http://blog.obeautifulcode.com/R/How-R-Searches-And-Finds-Stuff/
