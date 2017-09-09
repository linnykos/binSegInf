## Synopsis: See if manual rejection sampling works for selected model testing.
n=20
numIntervals= 10
library(genlassoinf)
library(binSegInf)
library(MASS)
source("../main/justin/sim-helper.R")

## Try wild binary segmentation fused lasso (1 step)
## nsim = 1000
## nsamp = 1000
nsim=1000
nsamp=300
addsigma=.1
## method = "wbs"
## type="segment"
method = "wbs"
type="fixed"
cat(fill=TRUE)
pvs = mclapply(1:nsim, function(isim){cat('\r', isim, "out of", nsim)
    z1list = list() ## Original data / model / contrast
    mu = rep(0,n)
    y0 = mu + rnorm(n,0,1)

    ## ## WBS
    ## I0 = binSegInf::generate_intervals(length(y0), numIntervals)
    ## original.model = wildBinSeg_fixedSteps(y0,intervals=I0, numSteps=1)
    ## original.cp = original.model$cp * original.model$cp.sign

    ## Fused lasso
    noise0 = rnorm(n,0,0.1)
    original.model = dualpathSvd2(y=y0 + noise0,
                                  D=makeDmat(n,ord=1,type='tf'),
                                  maxsteps=1)
    original.cp = original.model$cp * original.model$cp.sign

    contrasts <- make_all_segment_contrasts(original.model)
    v = contrasts[[as.character(original.cp[1])]]
    v = v/sqrt(sum(v*v))

    Proj = cbind(v)%*%rbind(v)
    Proj.perp = diag(rep(1,n)) - cbind(v)%*%rbind(v)

    ## Generate new data
    z0 = MASS::mvrnorm(nsamp, rep(0,n), Proj%*%t(Proj))
    z1mat = apply(z0, 1, function(myrow) Proj.perp%*%(y0+noise0) + myrow) ## For Fused lasso
    z1list = sanity.check(lapply(1:nsamp, function(isim)z1mat[,isim]), Proj.perp, 1E-6)

    ## Sample intervals (or noises)
    Ilist = lapply(1:nsamp, function(isim)binSegInf::generate_intervals(n, numIntervals))

    ## ## Get WBS rejection sampled z's
    ## cond.z1list = Map(function(z1,I){
    ##     if(!any(I$starts <= abs(original.cp) & abs(original.cp)<I$ends)){
    ##         return(NULL) ## This line causes some zeros?
    ##     }
    ##     new.model = wildBinSeg_fixedSteps(z1,intervals=I,numSteps=1)
    ##     new.cp = new.model$cp * new.model$cp.sign
    ##     if(all.equal(original.cp,new.cp)==TRUE) return(z1) else return(NULL)
    ## },z1list,Ilist)

    ## Get Fused Lasso rejection sampled z's
    addednoiselist = lapply(1:nsim,function(isim){rnorm(n,0,0.1)})
    cond.z1list = Map(function(z1,addednoise){
        new.model = dualpathSvd2(z1+addednoise,maxsteps=1,D=makeDmat(n,ord=1,type="tf"))
        new.cp = new.model$cp * new.model$cp.sign
        if(all.equal(original.cp,new.cp)==TRUE) return(z1) else return(NULL)
    },z1list,addednoiselist)

    ## Form the p-value by taking the quantile
    cond.z1list = cond.z1list[which(!sapply(cond.z1list,is.null))]
    nsurvive = length(cond.z1list)
    a = lapply(cond.z1list, function(z1) sum(v*z1))
    pv = sum(a>sum(v*y0))/nsurvive

    return(pv)

},mc.cores=3)

qqunif(unlist(pvs))

## pdf("~/Desktop/example.pdf",width=5,height=5)
## qqunif(c(unlist(pvs)))
## graphics.off()

