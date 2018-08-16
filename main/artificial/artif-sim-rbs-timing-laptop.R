## Synopsis: measuring time of randomized Binary segmentation inference. What
## factors matter? bits? number of steps? stopping rule? What part takes so
## long?


## Task 1 : Effect of min.num.things on stability of p-values (make sure same
## data is used for everything) <- Didn't show that much of a difference. What
## about for n=2000? Is there an obvious reason for the variance to change?
source("../main/artificial/artif-helpers.R")
minlist = c(10,20,30,40)
n = 200
sigma = 1
lev = 3
mn = c(rep(0,n/2), rep(lev,n/2))
y = mn + rnorm(n, 0, sigma)
bits = 1000
max.numIS = 2000
numIntervals = length(y)
min.num.things = 30
mc.cores = 8
nsim = 100
output.list = list()
for(imin in 1:length(minlist)){
    min.num.things=minlist[imin]
    printprogress(min.num.things, minlist, "minimum number of valid importance sampling replicates")
    cat(fill=TRUE)
    outputs = mclapply(1:nsim,function(isim){
        printprogress(isim,nsim)
        do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma,
                               postprocess=TRUE, locs=1:length(y), numIS=100,
                               sigma.add = 0.2, bits=bits,
                               inference.type="pre-multiply",
                               max.numIS=max.numIS,
                               min.num.things=min.num.things, verbose=FALSE)
    }, mc.cores=4)
    output.list[[imin]] = outputs
    cat(fill=TRUE)
}


## Task 2: min.num.things=10 stable enough? (making sure same data and same
## additive noise is used for everything)
source("../main/artificial/artif-helpers.R")
outputdir="../output"
n = 2000
sigma = 1
lev = 3
mn = rep(0,n)
mn[100:200] = mn[1100:1200] = mn[1800:2000] = lev
set.seed(99)
y = mn + rnorm(n, 0, sigma)
bits = 1000
max.numIS = 2000
numIntervals = length(y)
mc.cores = 8
output.list = list()
min.num.things=10
cat(fill=TRUE)
nsim = 20
sigma.add = 0.2
set.seed(2231)
added.noise = rnorm(n, 0, sigma.add)
outputs = mclapply(1:nsim,function(isim){
    printprogress(isim,nsim)
    do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma,
                     postprocess=TRUE, locs=1:length(y), numIS=100,
                     sigma.add = sigma.add, bits=bits,
                     inference.type="pre-multiply",
                     max.numIS=max.numIS,
                     min.num.things=min.num.things, verbose=FALSE,
                     added.noise = added.noise)
}, mc.cores=4)
cat(fill=TRUE)
save(outputs, file=file.path(outputdir, "rbs-stability.Rdata"))


## Load the output and analyze
load(file=file.path(outputdir, "rbs-stability.Rdata"))
hist(sapply(output.list[[1]],function(a)a$pvs)[1,], xlim=c(0,1))
hist(sapply(output.list[[2]],function(a)a$pvs)[1,])
sapply(output.list[[3]],function(a)a$pvs)
sapply(output.list[[3]],function(a)a$locs.all)


## Task 3: What takes the most time in a single n=2000 rbs simulation?
source("../main/artificial/artif-helpers.R")
outputdir="../outpt"
n = 2000
sigma = 1
lev = 5
mn = rep(0,n)
mn[100:200] = mn[1100:1200] = mn[1800:2000] = lev
## mn[10:20] = mn[110:120] = mn[180:200] = lev
set.seed(99)
y = mn + rnorm(n, 0, sigma)
bits = 2000
max.numIS = 2000
numIntervals = length(y)
mc.cores = 8
output.list = list()
min.num.things=10
cat(fill=TRUE)
nsim = 20
sigma.add = lev/3
set.seed(2231)
added.noise = rnorm(n, 0, sigma.add)
numIS=1
source("../main/artificial/artif-helpers.R")
Rprof(tmp <- tempfile())
output = do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma,
                          postprocess=TRUE, locs=1:length(y), numIS=numIS,
                          sigma.add = sigma.add, bits=bits,
                          inference.type="pre-multiply", max.numIS=max.numIS,
                          min.num.things=min.num.things, verbose=TRUE,
                          added.noise = added.noise)
Rprof()
summaryRprof(tmp)

output





## Task 3 result 1: Did I improve timing at all? from 11.5 to 6.9 by about 40%.
timing = do_rbs_inference()
## Task 3 result 2: Does this improve timing for n=2000? from OOO to OOO by about OO%.


## Task: run time by varying a few things.
## (1) by min.num.things [not done]
## (1-1) Is min.num.things=10 stable enough? How stable? [running]
## (2) by max.num.IS [not done]
## (3) by number of steps (under all else equal) [not done]
## (4) with/without stopping rule [not done]
## (5) Dependence on sample size?

