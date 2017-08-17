## Synopsis: See if manual sampling works.

## Setting:
n=10
numIntervals= 10
library(genlassoinf)
library(binSegInf)
library(MASS)


## Helpers
do.wildbinseg <- function(y){
    I = binSegInf::generate_intervals(length(y), numIntervals)
    original.model = binSegInf::wildBinSeg_fixedSteps(y,intervals=I, numSteps=1)
}
do.fl <- function(y){
    original.model = genlassoinf::dualpathSvd2(y,
                                               D=genlassoinf::makeDmat(length(y),ord=1,type='tf'),
                                               maxsteps=1)
}

## Try wild binary segmentation fused lasso (1 step)
nsim = 1000
nsamp = 200
addsigma=.1
## method = "wbs"
## type="segment"
method = "wbs"
type="fixed"
cat(fill=TRUE)
pvs = mclapply(1:nsim, function(isim){
    cat('\r', isim, "out of", nsim)
    z1list = list()

    ## Original data / model / contrast
    mu = rep(0,n)
    y0 = mu + rnorm(n,0,1)
    if(method=="wbs"){
        original.model = do.wildbinseg(y0)
    }
    if(method=="fl"){
        y1 = y0 + rnorm(n,0,addsigma)
        original.model = do.fl(y0)
    }
    original.cp = original.model$cp * original.model$cp.sign

    if(type=="fixed"){
        v = runif(n)
        v = v/sqrt(sum(v*v))
    }
    if(type=="segment"){
        contrasts <- make_all_segment_contrasts(original.model)
        v = contrasts[[as.character(original.cp[1])]]
    }
    Proj = cbind(v)%*%rbind(v)
    Proj.perp = diag(rep(1,n)) - cbind(v)%*%rbind(v)

    ## Generate new data
    z0 = MASS::mvrnorm(nsamp, Proj%*%mu, Proj%*%t(Proj))
    z1mat = apply(z0, 1, function(myrow) Proj.perp%*%y0 + myrow)
    z1list = lapply(1:nsamp, function(isim)z1mat[,isim])

    ## Rejection sample
    if(method=="wbs"){

        ## Sample intervals
        Ilist = lapply(1:nsamp, function(isim)binSegInf::generate_intervals(n, numIntervals))
        ## Get rejection sampled z's
        cond.z1list = Map(function(z1,I){
            new.model = binSegInf::wildBinSeg_fixedSteps(z1,intervals=I,numSteps=1)
            new.cp = new.model$cp * new.model$cp.sign
            if(all.equal(original.cp,new.cp)==TRUE) return(z1) else return(NULL)
        },z1list,Ilist)
    }
    if(method=="fl"){
        ## Get rejection sampled z's
        cond.z1list = lapply(z1list, function(z1){
            y1 = y0 + rnorm(n,0,addsigma)
            new.cp = do.fl(y1)
            if(all.equal(original.cp,new.cp)==TRUE) return(z1) else return(NULL)
        })
    }

    ## Form the quantile
    cond.z1list = cond.z1list[which(!sapply(cond.z1list,is.null))]
    nsurvive = length(cond.z1list)
    a = lapply(cond.z1list, function(z1) sum(v*z1))
    pv = sum(a>sum(v*y0))/nsurvive

    return(pv)

},mc.cores=2)

## qqunif(unlist(pvs))


