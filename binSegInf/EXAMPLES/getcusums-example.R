n = 30
sigma = 1
lev1 = 0
lev2 = 5
mn = c(rep(lev1,n/3), rep(lev2,n/3), rep(lev1,n/3))
set.seed(0)
y = mn + rnorm(n,0,sigma)
getcusums(1,n,y)
getcusums(1,2*n/3,y)
