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
obj = wildBinSeg_fixedThresh(y,thresh, numIntervals)
if(length(obj$cp)==0)return(NULL)

## Collect a polyhedron
poly <- polyhedra(obj)

## Make a test contrast
v = make_all_segment_contrasts(obj)[[1]]

## Get p-value conditioned on drawn interval set
pv1 = pvalue(y, poly, v, sigma)

## Get /randomized/ p-value /not/ conditional on interval set
pv2 <- randomized_wildBinSeg_pv(y = y, v = v,thresh=thresh,
                                numIntervals = numIntervals,
                                nsim.is=10, sigma=sigma)
