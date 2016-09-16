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
# Replacement of all usages of pval.fl1d to poly.pval; 
pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
}

## Setting 1
## n = 4
## set.seed(0)
## y = c(rep(0,n/2),rep(5,n/2)) + rnorm(n,0,1)
## thresh = .1

## Setting 2
## n = 60
## sd = .5
## mn = c(rep(-3, n/4), rep(2, n/4), rep(-1, n/4), rep(1, n/4))
## set.seed(1)
## y = mn + rnorm(n,0,sd)

## Setting 3
n = 20 
sigma=1
lev1=0
lev2=3
thresh = 2
nsim=1000

pvals.list = list()
pvals = matrix(NA,nrow=nsim,ncol=n)
for(isim in 1:nsim){
    cat('\r', isim, "out of", nsim)
    ## Set it up
    y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
    slist = elist = blist = Blist = zlist = Zlist =
        matrix(NA, nrow = n, ncol = 2^8)

    ## Do Binary Segmentation
    binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n)
    
    ## Bundle into one thing
    bs.output = list(slist = slist, elist = elist, blist = blist, Blist = Blist,
                     zlist = zlist, Zlist = Zlist, y = y, thresh = thresh)
    
    ## Get Gamma matrix and u vector
    Gu = get.polyhedron(bs.output, verbose=F)
    
    ## Conduct segment TG tests
    test.b.list = sort(collapse(blist))
    pvals[isim, test.b.list] = sapply(test.b.list,
                                 function(test.b){return( pval.fl1d(y = y,
                                                                    G = Gu$G,
                                                                    dik = make.v(test.b, bs.output),
                                                                    sigma = sigma,
                                                                u = Gu$u))})
}
pvals[test.b,]
hist(pvals[,10])
unif.sample = seq(from=0,to=1,length=nsim)
qqplot(pvals[,10], unif.sample)

replicate(100, do.once())

## If I want R to find stuff well, read this:
## http://blog.obeautifulcode.com/R/How-R-Searches-And-Finds-Stuff/
