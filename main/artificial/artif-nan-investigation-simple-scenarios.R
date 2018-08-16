## Synopsis: Trying out a one-jump signal with a really high signal size

# Data directory
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))

pvslist = list()
nsim = 20
numIS = 100
levs = c(5,10,30,50,100,200)
lev=5
for(ilev in 1:length(levs)){
    print("level")
    lev = levs[ilev]

    lev=100
    print(lev)
    pvslist[[isim]] = sapply(1:nsim, function(isim){
        printprogress(isim, nsim)

        n = numIntervals = 20
        mn = c(rep(0,n/2), rep(lev,n/2))
        sigma = 1
        y = mn + rnorm(n, 0, sigma)
        max.numSteps = 15
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=max.numSteps,
                                  inference.type='none')
        cumsum.y = cumsum(y)

        ## Collect the IC information and polyhedron
        ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
        ic_poly = ic_obj$poly
        stoptime  = ic_obj$stoptime

        ## Do the tests
        cp = g$cp[1:stoptime]
        cp.sign = g$cp.sign[1:stoptime]
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y), scaletype = "unitnorm")
        iv = 1
        v = vlist[[iv]]
        cumsum.v = cumsum(v)
        for(bits in c(20,50,100,200)){
            bits=100000
            for(isim in 1:3){
                pv = randomize_wbsfs(v=v, winning.wbs.obj=g,
                             sigma=sigma,
                             numIS=numIS,
                             inference.type=inference.type,
                             cumsum.y=cumsum.y,
                             cumsum.v=cumsum.v,
                             stop.time=stoptime+consec,
                             ic.poly=ic_poly,
                             improve.nomass.problem=improve.nomass.problem,
                             bits=bits, mc.cores=4)
            }
        }
    })
}


## How much longer does poly.pval() take with increasing number of bits? Answer: not longer at all?
bitslist = c(50,100,200,400,1000,2000,3000, 5000, 10000)
times.list = list()
for(isim in 1:nsim){
    print("isim is")
    print(isim)
    n = numIntervals = 20
    mn = c(rep(0,n/2), rep(lev,n/2))
    sigma = 1
    y = mn + rnorm(n, 0, sigma)
    max.numSteps = 15
    times = sapply(bitslist, function(bits){
        print("bits is")
        print(bits)
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=max.numSteps)
        poly = polyhedra(obj=g$gamma, u=g$u)
        ## Calculate TG denom and numer directly
        a = microbenchmark({
            pv = poly.pval2(v=v,y=y,sigma=sigma, poly=poly, bits=bits)
        }, times=100)
        mytime = print(a)$median
    })
    names(times) = bitslist
    times.list[[isim]] = times
}

ynew = rep(0,n) + rnorm(n, 0, sigma)
pv = poly.pval2(v=v,y=ynew,sigma=sigma, poly=poly)



z = sum(v*y)
vv = sum(v^2)
sd = sigma*sqrt(vv)
G = poly$gamma
u = poly$u
Gv = G %*% v
Gv[which(abs(Gv)<1E-15)] = 0
rho = Gv / vv
vec = (u - G %*% y + rho*z) / rho

which(is.nan(vec))
rho[which(is.nan(vec))]
(u - G %*% y + rho*z)[which(is.nan(vec))]
## vec = vec[rho!=0]
hist(vec)
vec[rho>0]
vec[rho<0]

vlo = suppressWarnings(max(vec[rho>0]))
vup = suppressWarnings(min(vec[rho<0]))
pv = tnorm.surv(z,0,sd,vlo,vup,bits)

matplot(do.call(cbind,times.list), type='l')

## What is that number? for our example?


## Make a single simulation example where you will get location 1859 (or 1862;
## not sure why it would differ this time around)
set.seed(0)
sigma = sd(y.orig[1:200])
y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)])

numIntervals=length(y)
consec=2
postprocess=TRUE
better.segment=TRUE
locs=1862
numIS=100; inference.type="pre-multiply"
improve.nomass.problem=TRUE; bits=1000; write.time=TRUE
verbose=TRUE
how.close=5

max.numSteps = 20
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=max.numSteps,
                          inference.type='none')
cumsum.y = cumsum(y)

## Collect the IC information and polyhedron
ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
ic_poly = ic_obj$poly
stoptime  = ic_obj$stoptime

if(verbose) cat("stoptime is", stoptime, fill=TRUE)

## Check for flag
if(ic_obj$flag!="normal" ){
    return(NULL)
}

## Extract changepoints from stopped model and declutter
cp = g$cp[1:stoptime]
cp.sign = g$cp.sign[1:stoptime]

if(postprocess){
    cpobj = declutter_new(cp, cp.sign, how.close=how.close)
    cp = abs(cpobj)
    cp.sign = sign(cpobj)
}

## # ## Plot things
## plot(y)
## abline(v=cp, col='purple',lwd=2)
## abline(v=g$cp, col="grey80")
## ## abline(v=g$cp[1:ic_obj$stoptime], col='blue', lwd=2)
## text(x=g$cp+3, y=rep(1,length(g$cp)), label = 1:length(g$cp))


## Form contrasts
if(better.segment){
    vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g, cps=cp)
} else {
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y))
}

## Retain only the changepoints we want results from:
print((abs(as.numeric(names(vlist)))))
retain = which((abs(as.numeric(names(vlist))) %in% locs))
if(length(retain)==0) return(list(pvs=c(), null.true=c()))

## Temporarily added
if(length(retain)!=0) browser()

## Calculate the p-values
vlist = vlist[retain]
iv = 1
v = vlist[[iv]]
cumsum.v = cumsum(v)
pv = randomize_wbsfs(v=v, winning.wbs.obj=g,
                                      sigma=sigma,
                                      numIS=numIS,
                                      inference.type=inference.type,
                                      cumsum.y=cumsum.y,
                                      cumsum.v=cumsum.v,
                                      stop.time=stoptime+consec,
                                      ic.poly=ic_poly,
                                      improve.nomass.problem=improve.nomass.problem,
                                      bits=bits, mc.cores=4)

## Why is the pvalue NaN? What changes it?
