n=30
lev=3
mn = c(rep(0,n/3),rep(lev,n/3),rep(0,n/3))
y = mn + rnorm(n,0,1)
a = circBinSeg(y)

## 
