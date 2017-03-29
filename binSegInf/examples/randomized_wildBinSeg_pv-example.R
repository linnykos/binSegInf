## Example settings
lev=3
thresh=2
numIntervals=10
n=10
sigma=1

## Generate data
mn <- rep(c(0,lev), each=n/2)
set.seed(2)
y <- mn + rnorm(n, 0, sigma)

## Run method /once/, collect things
set.seed(3)
obj = wildBinSeg(y,thresh, numIntervals)
if(length(obj$cp)==0)return(NULL)

## Collect a polyhedron
poly <- polyhedra(obj)

## Make a test contrast
v = make_all_segment_contrasts(obj)[[1]]

## Get p-value conditioned on drawn interval set
pv1 = pvalue(y, poly, v, sigma)

## Get /randomized/ p-value /not/ conditional on interval set
yo <- function(i){
    print(i)
pv2 <- randomized_wildBinSeg_pv(y = y, v = v,
                                numIntervals = numIntervals,
                                nsim.importance.sampling=1000)
return(pv2)
}

pv2s = sapply(1:100, yo)
hist(unlist(pv2s))
print(pv1)
print(pv2)
