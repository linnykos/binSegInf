## Synopsis: one-jump inference examples
source("../main/wbs-tests/plot-helpers.R")
source("../main/wbs-tests/sim-helpers.R")
outputdir = "../output"

## ## One-jump, one-step simulations

## ns = c(10,20,50,100,200)
## runtimes2 = list()
## for(i.n in 1:length(ns)){
##     n = ns[i.n]
##     nsim=1
##     numSteps=1
##     print(n)
##     runtimes2[[i.n]] = microbenchmark({
##         dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump)
##     }, times=5)
## }

## microbenchmark({
##         n=10
##         dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump)
##     }, times=5)

## nsim=1
## n=50
## microbenchmark({
##     print(n);print(nsim)
##     dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump)
## }, times=5)


## nsim=1
## numSteps=4
## system.time({
## dosim(lev=0, n=n, nsim=nsim, numSteps=numSteps, randomized=TRUE, numIS=100, meanfun=onejump, mc.cores=1)
## })


## load("../output/runtimes.Rdata")
## runtimes
## (runtimes[[1]])$time

## pdf(file="~/Desktop/runtimes.pdf", width=5, height=5)
## format(as.numeric((runtimes[[1]])$time), units="Mb",standard="legacy")
## mediantimes = sapply(1:5,function(ii){(as.numeric((runtimes[[ii]])$time)/1000000000)[4] })
## plot(mediantimes~ns, type='o', ylab="seconds", xlab="n")
## graphics.off()


## From laptop
source("../main/wbs-tests/sim-helpers.R")
lev=0
n=60
meanfun=onejump
nsim=1
numIS=100
randomized=TRUE
numSteps=1


## Rprof("~/Desktop/timing.out")
## pv = dosim(lev, n, meanfun, nsim, numSteps, numIS=numIS,
##            randomized, mc.cores=1, numIntervals=n, inference.type = "pre-multiply")
## Rprof(NULL)
## a=summaryRprof("~/Desktop/timing.out")
## head(a$by.total,30)



## b4 = summaryRprof("~/Desktop.out")
## b$by.self
## b$by.total


## Time things, by sample size
ns = c(10,20,50,100,200,500,1000,2000)
## durs = list()
durs.after.change = list()
for(ii in 1:length(ns)){
    print(ii)
    n = ns[ii]
    durs.after.change[[ii]] = microbenchmark({
        pv = dosim(lev, n, meanfun, nsim, numSteps, numIS=numIS,
                   randomized, mc.cores=1, numIntervals=n, inference.type = "pre-multiply")
    },times=5)
    outputdir="../output"
    filename = file.path(outputdir, "after-change-durations.Rdata")
    save(durs.after.change, file=filename)
}
times = sapply(durs, function(a)(a$time)[3]/1000000000)
times.after = sapply(durs.after.change, function(a)(a$time)[3]/1000000000)



filename = "after-change-durations.Rdata"
filename = "before-change-durations.Rdata"
outputdir = "~/Desktop"


times = sapply(durs, function(a)(a$time)[3]/1000000000)
times.after = sapply(durs.after.change, function(a)(a$time)[3]/1000000000)

load(file.path(outputdir, filename))
w=10; h=5
pdf("~/Desktop/times.pdf",width=w, height=h)
ns = c(10,20,50,100,200,500,1000,2000)
ns = ns[1:7]
par(mfrow=c(1,2))
for(iplot in 1:2){
    if(iplot==1){
        plot(times~ns, type='o', lwd=2)
    } else {
        plot(times~ns, type='o', xlim=c(0,2100), ylim = c(0,9000), lwd=2)
    }
    points(times.after~ns, type='o', col='red')

    dat = data.frame(times=times.after, n1=ns, n2=ns^2)
    g = lm(times~., data=dat)
    summary(g)
    x = 1:2000
    fitted = predict(g, newdata=data.frame(n1=x, n2=x^2))
    lines(fitted~x, col='grey')
}
graphics.off()
## Make the plots
plot(times, type='o')
points(times.after, type='o', col='red')

filename = "before-change-durations.Rdata"
filename = "after-change-durations.Rdata"
load(file.path(outputdir, filename))



## Time things, by sample size
ns = c(10,20,50,100,200,500,1000,2000)
## durs = list()
durs.after.change = list()
for(ii in 1:length(ns)){G
    print(ii)
    n = ns[ii]
    durs.after.change[[ii]] = microbenchmark({
        pv = dosim(lev, n, meanfun, nsim, numSteps, numIS=numIS,
                   randomized, mc.cores=1, numIntervals=n, inference.type = "pre-multiply")
    },times=5)
}

times = sapply(durs, function(a)(a$time)[3]/1000000000)
times.after = sapply(durs.after.change, function(a)(a$time)[3]/1000000000)

## Make the plots
w=10; h=5
pdf("~/Desktop/times.pdf",width=w, height=h)
ns = ns[1:5]
par(mfrow=c(1,2))
for(iplot in 1:2){
    if(iplot==1){
        plot(times~ns, type='o', lwd=2)
    } else {
        plot(times~ns, type='o', xlim=c(0,2100), ylim = c(0,9000), lwd=2)
    }
    points(times.after~ns, type='o', col='red')

    dat = data.frame(times=times.after, n1=ns, n2=ns^2)
    g = lm(times~., data=dat)
    summary(g)
    x = 1:2000
    fitted = predict(g, newdata=data.frame(n1=x, n2=x^2))
    lines(fitted~x, col='grey')
}
graphics.off()


