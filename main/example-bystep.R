library(devtools)
load_all()

n=100
numsteps = n-1
nsim=100
sigma=.5
p = rep(NA,nsim)
for(isim in 1:nsim){
    y = rnorm(n,0,sigma)
    a = binseg.by.size(y, numsteps,verbose=TRUE)
    print(all(a$G%*%y>0))
}



v = make.v(a$B[1],a$B,a$Z,10)



