load_all()

## Example
n = 12
set.seed(0) 
y = c(rnorm(n/3,0,.5),rnorm(n/3,3,.5), rnorm(n/3,5,.5))
numsteps = n-1
plot(y)
a = binseg.by.size(y, numsteps, verbose=TRUE)

## Binary segmentation
numsteps = 2
nsim=100
par(mfrow=c(2,5))
sigma=.5
p = rep(NA,nsim)
for(isim in 1:nsim){
    y = rnorm(10,0,1)
    a = binseg.by.size(y, numsteps,verbose=FALSE)
    v = make.v(a$B[1],a$B,a$Z,10)
    print(G%*%y>0)

    ## p[isim] = pval.fl1d(y = y,
    ##                     G = a$G,
    ##                     dik = v,
    ##                     sigma = sigma,
    ##                     u = a$u)
}
