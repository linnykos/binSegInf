## Synopsis: Investigating various things with respect to WBS speed

## How long does n=200 take, in overhead? Compared to a segment test
library(microbenchmark)
type = match.arg(type)
if(is.null(visc)) visc=1:n
sigma = 1
n = 200
meanfun = fourjump
nsim = 2000
numIntervals = n
numSteps = 4
visc.firstjump = n/5 + c(-1,0,1)
levs = c(0.1,0.2,0.3,0.4,0.5, 0.75)[whichlev]
mn = meanfun(lev=lev,n=n)
y = mn + rnorm(n, 0, sigma)
cumsum.y = cumsum(y)
inference.type = "pre-multiply"
improve.nomass.problem = TRUE


    sigma = 1
    mn = meanfun(lev=lev,n=n)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    inference.type = "pre-multiply"
    improve.nomass.problem = TRUE

numIntervals=n
numSteps=4
numIS=100
bits=10000
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
    poly.wbs = polyhedra(obj=g$gamma, u=g$u)
    vlist <- make_all_segment_contrasts(g)
max.numIS = 2000

## Get the p-values
pvs = sapply(vlist, function(v){
    print(v)
    cumsum.v = cumsum(v)
    pv = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g, sigma=sigma,
                                          numIS=numIS, inference.type=inference.type,
                                          cumsum.y=cumsum.y,cumsum.v=cumsum.v,
                                          improve.nomass.problem=improve.nomass.problem,
                                          bits=bits, max.numIS=max.numIS
                                          ))
})

## End of this Synopsis



