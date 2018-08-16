## Synopsis: Do a very specific analysis of bits=1000,2000,3000 run times
## /exactly where it occurs/, at poly.pval(). See how the runtime for
## depends on n

## Generate data
source("../main/artificial/artif-helpers.R")
outputdir="../output"
n = 2000
## n = 200
sigma = 1
lev = 5
mn = rep(0,n)
mn[100:200] = mn[1100:1200] = mn[1800:2000] = lev
## mn[10:20] = mn[110:120] = mn[180:200] = lev
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
sigma.add = lev/3
set.seed(2231)
added.noise = rnorm(n, 0, sigma.add)
numIS=1
## Fit model and get IC information
n = length(y)
max.numSteps=10;consec=2
new.noise = rnorm(length(y),0,sigma.add)
h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=max.numSteps)
ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma+sigma.add, type="bic")
stoptime = ic_obj$stoptime

## Collect stopped model and postprocess
cp = h.fudged$cp[1:stoptime]
cp.sign = h.fudged$cp.sign[1:stoptime]
vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
## vlist = vlist[1]
vlist = vlist[2]
v=vlist[[1]]

## Append the IC polyhedron
ic.poly = ic_obj$poly
poly = polyhedra(h.fudged)
orig.fudged.obj = h.fudged
poly$gamma = rbind(poly$gamma, ic.poly$gamma) ## I don't like this rbind..
poly$u = c(poly$u, ic.poly$u)
numSteps=NA
premult = polyhedra.bsFs(orig.fudged.obj,
                         inference.type="pre-multiply",
                         new.noise=new.noise, v=v,
                         numSteps=numSteps, y=y)
## Append IC stopping to Gy, Gv, Gw
ic.Gy = as.numeric(ic.poly$gamma%*%y)
ic.Gv = as.numeric(ic.poly$gamma%*%v)
ic.Gw = as.numeric(ic.poly$gamma%*%new.noise)
premult$Gy = c(premult$Gy, ic.Gy)
premult$Gv = c(premult$Gv, ic.Gv)
premult$Gw = c(premult$Gw, ic.Gw)
premult$u = c(premult$u, ic.poly$u)

objs.list = list()
times = list()
## bits.list = seq(from=100,to=3500, by=300)
bits.list = seq(from=2000,to=3500, by=200)
for(ibits in 1:length(bits.list)){
    bits = bits.list[ibits]
    printprogress(bits,bits.list, "bits"); cat(fill=TRUE)
    jj = 1
    objs=list()
    times[[ibits]] = microbenchmark({
        printprogress(jj, 100)
        objs[[jj]] = poly_pval_from_inner_products(Gy=premult$Gy,
                                                  Gv=premult$Gv, v=v, y=y,
                                                  sigma=sigma,
                                                  u=premult$u - premult$Gw,
                                                  bits=bits)
        jj=jj+1
    }, times = 100)
    objs.list[[ibits]] = objs
}

## Plot times
mediantimes = sapply(times, function(mytime)summary(mytime)$median)
plot(y=mediantimes, x=bits.list, type='o')
lm(mediantimes~bits.list)

## After excessive testing, we conclude that for $n=2000$, bits=2000 is a
## reasonable choice. There is a sharp transition at $bits=3000$-ish for the
## first contrast. There is no such sharp transition for the second contrast.
## Also doesn't exist for n=200.





## Another task: Does bits=100, bits=500 vs bits=1000 vs bits=2000 vs bits=3000
## make a difference in number of valid simulations?
source("../main/artificial/artif-helpers.R")
outputdir="../output"
n = 2000
## n = 200
sigma = 1
lev = 5
mn = rep(0,n)
mn[100:200] = mn[1100:1200] = mn[1800:2000] = lev
## mn[10:20] = mn[110:120] = mn[180:200] = lev
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
sigma.add = lev/3
set.seed(2231)
added.noise = rnorm(n, 0, sigma.add)
numIS=30
source("../main/artificial/artif-helpers.R")
library(microbenchmark)
bits.list = c(100, 200,300, 500, 1000)
outputs.list = list()
times = list()
for(ii in 1:length(bits.list)){
    bits = bits.list[[ii]]
    printprogress(bits, bits.list, "bits")
    cat(fill=TRUE)
    jj = 1
    outputs = list()
    outputs[[jj]] = do_rbs_inference(y=y, max.numSteps=10, consec=2, sigma,
                                     postprocess=TRUE, locs=1:length(y), numIS=numIS,
                                     sigma.add = sigma.add, bits=bits,
                                     inference.type="pre-multiply", max.numIS=max.numIS,
                                     min.num.things=min.num.things, verbose=TRUE,
                                     added.noise = added.noise, mc.cores=1,
                                     improve.nomass.problem=FALSE,
                                    return.more.things=TRUE)
    outputs.list[[ii]] = outputs
}
