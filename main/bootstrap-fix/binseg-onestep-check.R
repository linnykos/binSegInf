library(binSegInf)
library(smoothmest)
source("helpers.R")

##' Re-written one-step binary segmentation function
binseg_onestep <- function(y){
    
    # Polyhedron and
    G = matrix(0,2*n-4,n)
    M = B = matrix(0,n-1,n)

    for (i in 1:(n-1)) {
      M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
      B[i,] = sqrt(i*(n-i)/n) * M[i,]
    }
    
    # Binary segmentation
    j2 = which.max(abs(B %*% y))
    s2 = sign(sum(B[j2,]*y))
    G[1:(n-2),] = t(s2*B[j2,] - t(B[-j2,]))
    G[(n-1):(2*n-4),] = t(s2*B[j2,] + t(B[-j2,]))
    return(list(gamma=G, u=rep(0, nrow(G)),
                y=y, cp=j2, cp.sign=s2))
}

## Compare output and polyhedra with existing code
nrep = 100
n = 100
for(irep in 1:nrep){## Generate data
    set.seed(irep)
    print(irep)
    y = rep(0, n) + rnorm(n, 0, 1)

    ## Fit two one-step algorithms and form contrast
    obj1 = binSeg_fixedSteps(y, numSteps=1)
    obj2 = binseg_onestep(y)
    stopifnot(obj1$cp * obj1$cp.sign == obj2$cp * obj2$cp.sign)
    i = obj1$cp
    v = obj1$cp.sign*c(rep(-1/i,i),rep(1/(n-i),n-i))
    v = v/sqrt(sum(v*v))

    ## Compare polyhedra
    G1 = polyhedra(obj1)$gamma
    G2 = obj2$gamma

    ## Compare p-values of resulting segment test
    p1 = pval(y, G1, v)
    p2 = pval(y, G2, v)

    ## Compare the polyhedra -- does it match for new data case by case?
    nsim=1000
    ynewmat = sapply(1:nsim, function(isim) y + rnorm(n,0,1))
    stopifnot(all(apply(ynewmat, 2, function(ynew){ sign(G1%*%ynew)==sign(G2%*%ynew)})))
    stopifnot(all.equal(p1, p2))
}

nrep=20000
out.gaus = sim(err.fun = function(n) rnorm(n,0,1), nrep=nrep)
out.unif = sim(err.fun = function(n) runif(n, min=-sqrt(12)/2, max=sqrt(12)/2), nrep=nrep)
out.lapl = sim(err.fun = function(n) rexp(n,rate=sqrt(2)) * 
            sample(c(-1,1),n,replace=TRUE), nrep=nrep)
outputdir = "../figures"
jpeg(file=file.path(outputdir, "fivethousand.jpg"), width=700, height=700)
qqunif(list(gaussian = out.gaus$p.bs, uniform = out.unif$p.bs, laplace=out.lapl$p.bs), cex=rep(.3, 3),
       col=c(1,2,3), legend=TRUE)
title(main="Laplace noise, 1-step BS")
graphics.off()

