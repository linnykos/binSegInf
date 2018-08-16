# Synopsis: See how long simulations are going to take.
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

##' Multiply the maximum to have noise*lev maximum
##' height.
outputdir = "../output"
list.of.cplists = list()
for(fac in 1:5){
    cplist = list()
    printprogress(fac, 5, "factor multiplied.")
    cat(fill=TRUE)
    for(seed in 1:100){
        printprogress(seed, 100)
        ## set.seed(seed)
        sigma = sd(y.orig[1:200]) * fac
        sigma.add = sigma*0.2
        ## y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)]) * fac
        y = newmn[-(1:200)] + bootstrap_sample

        ## Timing the inference.
        verbose=TRUE
        bits=1000
        numIS = 20
        min.num.things = 100
        how.close=5
        start.time = Sys.time()
        object = inference_bsFs(y=y, max.numSteps=15, consec=2,
                                sigma=sigma, postprocess= TRUE,
                                locs=1:length(y), numIS= numIS,
                                min.num.things=min.num.things,
                                inference.type="pre-multiply",
                                bits=bits, sigma.add=sigma.add,
                                verbose=verbose, start.time=start.time,
                                how.close=how.close,
                                return.more.things=TRUE)
        cplist[[seed]] = object
    }
    list.of.cplists[[fac]] = cplist

    ## See what the 1-length models are
    ## a = cplist[lapply(cplist, length)==1]
    ## a = cplist[lapply(cplist, length)==2]
    ## plot(a)

    ## Temporary load
    outputdir = "../output"
    load(filename=file.path(outputdir,"rbs-models.Rdata"))
    fac = 2
    cplist = list.of.cplists[[fac]]

    ## Make plot of results
    filename = paste0("artif-rbs-fac-",fac,".jpg") ## pdf
    ## pdf(file=file.path(outputdir, filename), width=7, height=7)
    jpeg(file=file.path(outputdir, filename), width=700, height=700)
    par(mfrow=c(2,1))

    ## Make barplots of stop times
    didntstop = sapply(cplist,function(a)all(is.na(a)))
    sum(didntstop)/100
    mytable = table(sapply(cplist[-didntstop],length))
    mytable = c(sum(didntstop), mytable)
    names(mytable)[1]="didn't-stop"
    barplot(mytable, main = "Model sizes stopped by 2-rise BIC")

    ## Clean some results
    cps = unlist(cplist)
    cps.sorted = sort(table(cps),decreasing=TRUE)
    cps.all = unlist(cps)
    cps.all = cps.all[which(!is.na(cps.all))]
    a = density(cps.all,width=20)

    ## Make plot
    plot(y, pch=16, col='grey50', ylim = c(-1.5,1.5))
    lines(newmn[-(1:200)], col='red', lwd=3)
    lines(a$y/max(a$y)*1.5~a$x, lwd=3)

    graphics.off()
}
save(list.of.cplists, file=file.path(outputdir, "rbs-models.Rdata"))




## Synopsis: Single example to see the number of importance sampling replicates required?
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
library(microbenchmark)

##' Multiply the maximum to have noise*lev maximum
##' height.
outputdir = "../output"
fac = 4
## set.seed(seed)
sigma = sd(y.orig[1:200]) * fac
sigma.add = sigma*0.2
resid.cleanmn2 = resid.cleanmn
newmn2 = newmn
## newmn2[1859:2109] = 2*newmn[1859:2109]## temporary, to boost the hardness of the
##                                       ## example? This failed; it isn't hard


set.seed(1)
y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn2[-(1:200)]) * fac
plot(y)
numIS = 10
min.num.things = 20
how.close = 5
bits=3000
start.time = Sys.time()
object = inference_bsFs(y=y, max.numSteps=15, consec=2,
                        sigma=sigma, postprocess= TRUE,
                        locs=1:length(y), numIS= numIS,
                        min.num.things=min.num.things,
                        inference.type="pre-multiply",
                        bits=bits, sigma.add=sigma.add,
                        verbose=verbose, start.time=start.time,
                        how.close=how.close,
                        return.more.things=TRUE,
                        ## myloc=1859+c((-2):2)
                        myloc=923+c((-2):2)
                        )

## ## What does the v look like?
## v=object[[1]]
## plot(y)
## lines(v*30, lwd=3)




## How often does it find the first four jumps at fac=4?
len = length(cp)
mytime = list()
pvs.rbs = c()
objects = list()

for(ii in 1:5){
    set.seed(ii)
    start.time = Sys.time()
    objects[[ii]] = inference_bsFs(y=y, max.numSteps=10, consec=2,
                                   sigma=sigma, postprocess= TRUE,
                                   locs=1:length(y), numIS= numIS,
                                   min.num.things=min.num.things,
                                   inference.type="pre-multiply",
                                   bits=bits, sigma.add=sigma.add,
                                   verbose=verbose, start.time=start.time,
                                   how.close=how.close,
                                   return.more.things=TRUE,
                                   whichv=1)
    cat(fill=TRUE)
}



