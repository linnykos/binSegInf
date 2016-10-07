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

n = 20
sigma=1
lev1=0
lev2=3
thresh = 2
nsim=1000
lambda=1
lev2s = seq(from=0, to=5, by=.5)

## Generate p-values
pvals.list = list()
for(ilev2 in 1:length(lev2s)){
    cat(ilev2, "out of", length(lev2s), 'level2s', '\n')
    pvals = matrix(NA, nrow=nsim, ncol=n)
    lev2 = lev2s[ilev2]
    for(isim in 1:nsim){
        cat('\r', isim, "out of", nsim, 'simulations')
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
        Gu = get.polyhedron(binseg.results = bs.output, lambda,verbose=F)
    
        ## Conduct segment TG tests
        test.b.list = sort(collapse(blist))
        pvals[isim, test.b.list] = as.numeric(sapply(test.b.list,
                                          function(test.b){return( pval.fl1d(y = y,
                                                                             G = Gu$G,
                                                                             dik = make.v(test.b, bs.output),
                                                                             sigma = sigma,
                                                                             u = Gu$u))}))
    }
    cat(fill=T)
    pvals.list[[ilev2]] = pvals
}

sim.settings = list(lev1=lev1,lev2=lev2,thresh=thresh,nsim=nsim,sigma=sigma, n=n)
save(pvals.list, sim.settings, file='onejump-by-lev2.Rdata')

## Plotting code


## Plotting code
load('onejump-by-lev2.Rdata')
pdf("~/Desktop/onejump-by-lev2.pdf")
par(mfrow=c(3,4))
lev2 = 3
set.seed(0)
y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
plot(y, pch = 16, main = "data and signal (example)")
lines(c(rep(lev1,n/2),rep(lev2,n/2)), col='red',lwd=2)
for(ilev2 in 1:length(lev2s)){
    hist(pvals.list[[ilev2]][,10],main="")
    title(main=bquote(delta == .(lev2s[ilev2])))
}
graphics.off()

## unif.sample = seq(from=0,to=1,length=nsim)
## qqplot(pvals[,10], unif.sample)

## replicate(100, do.once())

## If I want R to use environments to find stuff, read this:
## http://blog.obeautifulcode.com/R/How-R-Searches-And-Finds-Stuff/
