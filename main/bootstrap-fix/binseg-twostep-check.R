library(binSegInf)
library(smoothmest)
source("helpers.R")
library(genlassoinf)

n = 10
mn = rep(0, n)
set.seed(0)
y = mn + rnorm(n,0,1)

##' Re-written two-step binary segmentation function
binseg_twostep <- function(y){
    
    # Polyhedron and
    G1 = matrix(0,2*n-4,n)
    G2 = matrix(0,2*n-6,n)
    M = B = matrix(0,n-1,n)

    for (i in 1:(n-1)) {
      M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
      B[i,] = sqrt(i*(n-i)/n) * M[i,]
    }
    
    # First step binary segmentation
    cp1 = which.max(abs(B %*% y))
    s1 = sign(sum(B[cp1,]*y))
    G1[1:(n-2),] = t(s1*B[cp1,] - t(B[-cp1,]))
    G1[(n-1):(2*n-4),] = t(s1*B[cp1,] + t(B[-cp1,]))

    ## Manually form the second step binary segmentation
    D = genlassoinf::makeDmat(m=n, order=0)
    D1 = D[-cp1,,drop=FALSE]
    alin = MASS::ginv(t(D1))
    base.coef = do.call(rbind,
                        lapply(1:nrow(alin), function(ii){
                            sbe = get_sbe(alin[ii,])
                            v = make.base.contrast.from.sbe(sbe[["s"]],sbe[["b"]],sbe[["e"]],n)
                            return(v)
                        }))
    weights.binseg = apply(base.coef, 1, function(myrow)1/sqrt(sum(myrow^2)))
    binseg.coef = base.coef * weights.binseg
    B = binseg.coef
    cp2 = which.max(abs(B %*% y))
    s2 = sign(sum(B[cp2,]*y))
    G2[1:(n-3),] = t(s2*B[cp2,] - t(B[-cp2,]))
    G2[(n-2):(2*n-6),] = t(s2*B[cp2,] + t(B[-cp2,]))
    gamma = rbind(G1,G2)

    ## Form the two contrasts
    vlist = make_all_segment_contrasts_from_cp(c(cp1,cp2), c(s1,s2), n, "unitnorm")
    v1 = vlist[[1]]
    v2 = vlist[[2]]

    ## Calculate p-values
    p1 = pval(y, gamma, v1)$p
    p2 = pval(y, gamma, v2)$p
    pvs = c(p1,p2)
    names(pvs) = c(cp1,cp2)
    
    return(pvs)
}

## out.new = binseg_twostep(y)
## out = binSeg_fixedSteps(y, numSteps=2)


## ## Two-step polyhedron is the same, entry-wise! 
## gamma.new = out.new$gamma
## gamma = polyhedra(out)$gamma
## all.equal(gamma.new, gamma)


## The two-step inference is as skewed as I originally found? Just checking one
## more time!
nrep = 20000
n = 100

start.time=Sys.time()
mn = rep(0, n)
pvslist.gaus = mclapply(1:nrep,function(irep){
    printprogress(irep,nrep,start.time=start.time)
    y = mn + rnorm(n,0,1)
    return(binseg_twostep(y))
},mc.cores=3)

pvslist.lapl = mclapply(1:nrep,function(irep){
    printprogress(irep,nrep,start.time=start.time)
    y = mn + rexp(n,rate=sqrt(2)) * sample(c(-1,1),n,replace=TRUE)
    return(binseg_twostep(y))
},mc.cores=3)

pvslist.unif = mclapply(1:nrep,function(irep){
    printprogress(irep,nrep,start.time=start.time)
    y = mn + runif(n, min=-sqrt(12)/2, max=sqrt(12)/2)
    return(binseg_twostep(y))
},mc.cores=3)

outputdir = "../output"
filename = "twostep-binseg.Rdata"
save(pvslist.gaus, pvslist.lapl, pvslist.unif,
     file=file.path(outputdir, filename))

outputdir="."
load( file=file.path(outputdir, filename))

outputdir="../figures"
jpeg(file=file.path(outputdir, "twostep-binseg.jpg"), width=700, height=700)
qqunif(list(gaussian = do.call(rbind, pvslist.gaus),
            uniform = do.call(rbind, pvslist.unif),
            laplace = do.call(rbind, pvslist.lapl)), legend=TRUE, cols=c(1:3))
title(main="Two-step binary segmentation")
graphics.off()
