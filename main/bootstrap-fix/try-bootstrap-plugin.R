## Synopsis: Using the bootstrap plugin of 5.2 in asympinf paper (but with our
## /known/ residuals instead of $Y-\bar Y$)
library(binSegInf)
library(smoothmest)
source("helpers.R")
outputdir = "../figures"
la("~/repos/binSegInf/binSegInf")
la("~/repos/genlassoinf/genlassoinf/")

## Run simulation
lapl <- function(n){ rexp(n,rate=sqrt(2)) * sample(c(-1,1),n,replace=TRUE)}
nrep = 3000##500000
out.lapl = simnew(n=100, nrep=nrep, err.fun=lapl, seed=NULL)
out.gaus = simnew(n=100, nrep=nrep, err.fun=rnorm, seed=NULL)

## Make plots
makejpg(outputdir, 
qqunif(out.lapl[11:14], cols=1:4)
qqunif(out.gaus[11:14], cols=1:4)

## How about null cases in nonzero jump signals?
rep = 100
n = 10
p.tg = p.plugin = list()
lev = 3
mn = c(rep(0, n/2), rep(lev, n/2))
nboot = 5000
for(irep in 1:nrep){## Generate data

    printprogress(irep, nrep)
    y = mn + rnorm(n, 0, 1)

    ## Fit two one-step algorithms and form contrast
    obj = binSeg_fixedSteps(y, numSteps=3)
    vlist = make_all_segment_contrasts(obj)
    tol = 1E-10

    ## Return only null contrasts
    ## whichnull = which(sapply(vlist, function(v) abs(sum(v*mn))<tol ))
    ## vlist = vlist[whichnull]
    ## if(length(vlist)==0)next

    ## Return only nonnull contrasts
    whichnonnull = which(!sapply(vlist, function(v) abs(sum(v*mn))<tol ))
    vlist = vlist[whichnonnull]
    if(length(vlist)==0)next

    ## Compare polyhedra
    G = polyhedra(obj)$gamma

    ## Compare p-values of resulting segment test
    p.tg[[irep]] = sapply(vlist, function(v){
        ## pval(y, G, v)$p
        poly.pval(y=y, G=G, u=rep(0,nrow(G)), v, sigma=1)$p
    })
    p.plugin[[irep]] = sapply(vlist, function(v)pval_plugin_wrapper(y, G, v, nboot=nboot))
    
    ## source('helpers.R')
    ## pval_plugin_wrapper(y, G, v, nboot=nboot)
}


p.plugin.5000 <- p.plugin
p.plugin.2000 <- p.plugin

p.plugin.5000
p.plugin.2000


p.plugin ## There are quite/ a few NaNs. What is going on there?
## p.plugins = (unlist(p.plugin))
p.plugin[[6]]<0.05/2
rejects = sapply(p.plugin, function(pvs){pvs<0.05/length(pvs)})
sum(unlist(rejects))/length(unlist(rejects))
(rejects)

## Why are NaNs happening?


## What about power?
source("helpers.R")
qqunif(list(plugin=unlist(p.plugin), tg=unlist(p.tg)), cols=c(1,2))


## There are a fair number of 1's. Why is this?

## Next up: do power comparison in 1,2,3 steps.
## Then next up: try to think abuot how this will go in the paper 


    ## pval_plugin = function(Vlo, Vup, vty, v, y) {



## Step 1: Store all the Vup/Vlo/vty/v

## Step 2: For each triplet Calculate the BSTG (bootstrap-TG) distribution.

## Currently, the bootstrap distribution is wrong. Trying to figure out why.

