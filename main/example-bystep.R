library(devtools)
load_all()

## Example
n = 12
set.seed(0) 
y = c(rnorm(n/3,0,.5),rnorm(n/3,3,.5), rnorm(n/3,5,.5))
numsteps = n-1
a = binseg.by.size(y, numsteps, verbose=TRUE)
binseg(y,0,verbose=TRUE)


## Binary segmentation
n=100
numsteps = n-1
nsim=100
sigma=.5
p = rep(NA,nsim)
for(isim in 1:nsim){
    y = rnorm(n,0,sigma)
    a = binseg.by.size(y, numsteps,verbose=FALSE)
    print(all(a$G%*%y>0))
    v = make.v(a$B[1],a$B,a$Z,10)
}


## Sanity check function
check.binseg = function(G,u)
    y = rnorm(10,0,1)
numsteps = 2
nsim=100
## par(mfrow=c(2,5))
sigma=.5



