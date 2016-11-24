library(devtools)
load_all()

n=12
numsteps = n-1
nsim=100
sigma=.5
y = rnorm(n,0,sigma)
a = binseg.by.thresh(y, 2,verbose=TRUE)
print(all(a$G%*%y>0))
