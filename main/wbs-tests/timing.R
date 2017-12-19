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





    getmedtime <- function(mytime){(mytime$time[3])/1000000000}
    medtimes= sapply(results, getmedtime)
