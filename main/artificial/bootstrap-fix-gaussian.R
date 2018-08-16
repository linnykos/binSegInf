
library(rmutil)

onesim_gaussian_replacement_rbs <- function(bits=1000,  verbose=FALSE, seed=NULL,
                                 max.numIS=2000, min.num.things=30){
    ## Generate Gaussian noise
    ## mn = fourjump(n=n, lev=lev)
    ## y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac
    if(!is.null(seed))set.seed(seed)
    ## subresid = subresd - mean(subresid) ## centering
    y = newmn + bootstrap_sample(resid.cleanmn)
    ## sigma = sd(resid.cleanmn)
    ## ## sigma.add=sigma*0.2
    ## sigma.add=0

    n=20
    lev=1
    mn = fourjump(n=n, lev=lev)
    resids2 <- rmutil::rlaplace(n,m=0,s=1/sqrt(2))
    ## y = mn + rnorm(n,0,1)
    y = mn + resids2
    sigma=1
    sigma.add=0

    ## Obtain inferences.
    postprocess = FALSE
    how.close = 1
    numIS = 10
    min.num.things = 30
    object = inference_bsFs(y=y, max.numSteps=10, consec=2, sigma=sigma,
                            postprocess= postprocess, locs=1:length(y), numIS= numIS,
                            min.num.things=min.num.things,
                            inference.type="pre-multiply", bits=bits,
                            sigma.add=sigma.add, verbose=verbose,
                            start.time=start.time, how.close=how.close,
                            max.numIS=max.numIS,
                            sim.options = list(retain.only.null.cases=TRUE,
                                               mn=mn,
                                               old=FALSE))
    return(object)
}

nsim=1000
    start.time = Sys.time()
    results = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim, start.time=start.time)
        onesim_gaussian_replacement_rbs(bits=5000,  verbose=FALSE,
                             max.numIS=2000, min.num.things=30)
    },mc.cores=4)

plotfilename="bootstrap-laplace.jpg"
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
qqunif(unlist(results))
graphics.off()


mc.cores = 4
nsim = 500
results.list = lapply(1:10, function(seed){
    printprogress(seed,100, "seed", fill=TRUE)
    start.time = Sys.time()
    results = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim, start.time=start.time)
        onesim_gaussian_replacement_rbs(bits=5000,  verbose=FALSE,
                             max.numIS=2000, min.num.things=30)
    })
    return(results)
})

filename = "bootstrap-gaussian-replacement.Rdata"
save(results, file=file.path(outputdir,filename))


filename = "bootstrap-gaussian-replacement.Rdata"
makeplotfilename <- function(filename, ext= ".jpg"){ gsub(".Rdata", ext, filename )}
plotfilename = makeplotfilename(filename)
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
par(mfrow=c(2,3))
for(ifac in 1:5){
    ## fac = c(1, 1.1, 1.2, 1.3)[ifac]
    result = results.list[[ifac]]
    qqunif(unlist(result))
    set.seed(ifac)
    subresid = resid.cleanmn[sample(length(resid.cleanmn), n)]
    ## plot(subresid, ylim=range(resid.cleanmn))
    ## abline(h=mean(subresid))
    title(main=round(sd(subresid),3))
}
graphics.off()



## Make qqplot
plotfilename = "bootstrap-resid-qqplot.jpg"
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
qqnorm(resid.cleanmn, main = "QQ plot of all residuals")
qqline(resid.cleanmn)
graphics.off()

## compare to laplace
plotfilename = "bootstrap-laplace-qqplot.jpg"
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
resids <- resid.cleanmn/sd(resid.cleanmn)
resids
resids2 <- rmutil::rlaplace(1000000,## length(resids)
                          m=0,s=1/sqrt(2))
sd(resids2)
qqnorm(resids2)
qqline(resids2)
graphics.off()


## Raw residual plot
lotfilename = "bootstrap-resid-first200.jpg"
s
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
qqnorm(resid.cleanmn[1:200], main = "QQ plot of all residuals")
qqline(resid.cleanmn[1:200])
graphics.off()

plotfilename = "bootstrap-resid-first400.jpg"
jpeg(file=file.path(outputdir, plotfilename), width=700,
     height=700)
qqnorm(resid.cleanmn[1:400], main = "QQ plot of all residuals")
qqline(resid.cleanmn[1:400])
graphics.off()
