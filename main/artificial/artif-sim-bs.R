## Synopsis: Same as artif-sim-rbs but trying plain binary segmentation

# Data directory
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
## Add bootstrapped residuals around a cleaned mean, with known sigma

## Simulation settings
onesim_bs <- function(y.orig, bits=1000, fac=1, verbose=FALSE){

    ## Add bootstrapped residuals around a cleaned mean, with known sigma
    ## fac=1
    mysd = sd(y.orig[1:200])*1.2
    ## plot(y.orig[1:200])
    ## qqnorm(y.orig[1:200])
    ## qqline(y.orig[1:200])
    ## bootstrap

    sigma = mysd * fac
    ## sigma.add = sigma*0.2
    sigma.add=0
    y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac

    how.close=5
    numIS=10
    min.num.things=30
    start.time = Sys.time()
    object = inference_bsFs(y=y, max.numSteps=15, consec=2,
                            sigma=sigma, postprocess= TRUE,
                            locs=1:length(y), numIS= numIS,
                            min.num.things=min.num.things,
                            inference.type="pre-multiply",
                            bits=bits, sigma.add=sigma.add,
                            verbose=verbose, start.time=start.time,
                            how.close=how.close,
                            ## myloc=1859+c((-2):2),
                            )
}

## Fix factor fac=1. See how to improve type I errors
start.time = Sys.time()
bits=2000
mc.cores=3
nsim=100
fac=1
results = mclapply(1:nsim, function(isim){
    printprogress(isim, nsim,
                  lapsetime = round(difftime(Sys.time(), start.time,
                                             units = "secs"), 2))
    ## Run a single result
    myresult = onesim_bs(y.orig, bits=bits, fac=fac, verbose=FALSE)
}, mc.cores=mc.cores)



qqunif(newmn[1"20"])

## This is the data segment
par(mfrow=c(1,2))
plot(y.orig[1:200])
qqnorm(y.orig[1:200])
qqline(y.orig[1:200])
graphics.off()

## This is the standard deviation's bootstrapped distribution.


