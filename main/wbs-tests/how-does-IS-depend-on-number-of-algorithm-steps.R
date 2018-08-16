## Synopsis: How does the total number of IS depend on the number of steps taken
## in the algorithm? (=amount of conditioning?)

outputdir = "../output"
source("../main/wbs-tests/sim-helpers.R")

## Generate data once
n = 200
lev = 1 ## lev=2
## lev=2
sigma = 1
mn = fourjump(lev=lev,n=n)
y = mn + rnorm(n, 0, sigma)
plot(y)
points(mn,lwd=3)
cumsum.y = cumsum(y)
inference.type = "pre-multiply"
improve.nomass.problem = TRUE
numIS = 100
max.numIS=2000
numIS=100
bits=3000
min.num.things = 30

## This is v[[1]]!!
v = rep(0,n)
v[1:120] = -1/120
v[121:n] = 1/(n-120)

set.seed(1)
intervals = intervals(n=n, numIntervals=n)

min.num.things=30
cumsum.v = cumsum(v)
cumsum.y = cumsum(y)

## Apply WBS
infolists = list()
maxNumSteps=6
for(numSteps in 1:maxNumSteps){
    printprogress(numSteps, maxNumSteps, "Steps taken in the algorithm")
    cat(fill=TRUE)
    g = wildBinSeg_fixedSteps(y, intervals=intervals, numSteps=numSteps)
    poly.wbs = polyhedra(obj=g$gamma, u=g$u)
    start.time = Sys.time()
    infolist = list()
    mc.cores=8
    nsim=100
    for(isim in 1:nsim){
        printprogress(isim,nsim, "total simulations",
                      lapsetime = round(difftime(Sys.time(), start.time,
                                                 units = "secs"), 2))
        infolist[[isim]] = randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                           numIS=numIS,
                                           inference.type=inference.type,
                                           cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                           improve.nomass.problem=improve.nomass.problem,
                                           bits=bits, max.numIS=max.numIS,
                                           return.more.things=TRUE, min.num.things=min.num.things,
                                           mc.cores=mc.cores, verbose=FALSE)
    }
    infolists[[numSteps]] = infolist
    ## Save results
    ## filename = "how-stable-is-rwbs.Rdata"
    filename = "how-does-IS-depend-on-conditioning.Rdata"
    save(list=ls(), file=file.path(outputdir, filename))
}

## Load
filename = "how-does-IS-depend-on-conditioning.Rdata"
load(file=file.path(outputdir, filename))

## How many replicates are required, for each case?
w=h=5;
## pdf(file=file.path(outputdir, "how-stable-is-rwbs-numdraws.pdf"), width=w, height=h)
pdf(file=file.path(outputdir, "how-does-IS-depend-on-conditioning.pdf"), width=w, height=h)
mean.draws = sd.draws = c()
for(ii in 1:maxNumSteps){
    infolist = infolists[[ii]]
    mean.draws[ii] = mean(sapply(infolist, function(a)a$numIS.cumulative))
    sd.draws[ii] = sd(sapply(infolist, function(a)a$numIS.cumulative))
}
mymat = cbind(mean.draws+sd.draws,mean.draws, mean.draws-sd.draws)
matplot(NA, ylim = range(as.numeric(mymat)), col='white',
        xlim = range(min.num.things.list), axes=FALSE,
        xlab = "Number of WBS algorithm steps taken",
        ylab = "total # of draws that were run")
axis(1);axis(2)
for(irow in 1:nrow(mymat)){
    myrow = mymat[irow,]
    print(myrow)
    lowtick = myrow[1]
    midtick =  myrow[2]
    uptick = myrow[3]
    x.to.plot = min.num.things.list[irow]
    segments(x0=x.to.plot,x1=x.to.plot,y0=lowtick,y1=uptick)
    points(x=x.to.plot,y=lowtick, pch="-")
    points(x=x.to.plot,y=midtick, pch=16)
    points(x=x.to.plot,y=uptick, pch="-")
}
title(main = "How many IS draws are needed for \n certain # of valid IS draws")
graphics.off()
