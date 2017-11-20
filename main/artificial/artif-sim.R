## Data directory
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

##' Multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1, n){
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}

## Simulation settings
onesim <- function(y.orig){

    n = length(y.orig)

    ## Generate nosie around the data
    sigma=.1
    y = newmn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)

    ## Do two types of inference: randwbs, randFL
    pvs.rwbs = do_rwbs_inference(y=y, max.numSteps=10, numIntervals=length(y),
                                 consec=2, sigma=sigma, postprocess=TRUE,
                                 better.segment=TRUE, locs=1:length(y),
                                 numIS=100, inference.type="pre-multiply",
                                 improve.nomass.problem=TRUE)
    rfl.time = microbenchmark({
        source(file=file.path("../main/artificial/artif-helpers.R"))
        pv.rfl = do_rfl_inference(y=y, max.numSteps=10,
                                  consec=2, sigma=sigma, postprocess=TRUE,
                                  locs=1:length(y), numIS=100, sigma.add = 0.05)
    }, times=1)
}
onesim(y.orig)


## One replicate ofRandomization is timed here
v = vlist[[1]]
            cumsum.v = cumsum(v)

rand.time = microbenchmark({
    numIS = 20
    pv.rand = randomize_wbsfs(v=v, winning.wbs.obj=g,
                    sigma=sigma,
                    numIS=numIS,
                    inference.type="pre-multiply",
                    cumsum.y=cumsum.y,
                    cumsum.v=cumsum.v,
                    stop.time=stoptime+consec,
                    ic.poly=NULL,
                    improve.nomass.problem=TRUE,
                    bits=1000)
})

fit.time = microbenchmark({
    g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                              inference.type='pre-multiply', cumsum.y=cumsum.y, cumsum.v=cumsum.v)
    poly_pval_from_inner_products(g$Gy,g$Gv, v, y, sigma, g$u, bits=1000)
},times=5)




## Time a single replicate of inference


levs=seq(from=0,to=3,length=n.levs)
nsims = seq(from=100,to=50,length=n.levs)
bootstrap=TRUE
reduce=TRUE
mc.scores = 1




## Plot things
pdf("~/Desktop/sample-data.pdf",width=10,height=10)
par(mfrow=c(2,2))
for(lev in 1:4){
mn <- coriell_mn
set.seed(0)
y <- mn(lev,n) + rnorm(n,0,std)
plot(y,ylim=c(-1,1), main = paste0("Signal size /stretched/ from snr=1 to ",lev),pch=16,cex=.5,col='grey50')
lines(mn(lev,n),col='red')
}
graphics.off()
