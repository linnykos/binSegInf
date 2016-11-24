require(wbs)

## Generate data
sigma = 1
lev1 = 0
lev2 = 3
n=100
f = c(rep(lev1,n/2), rep(lev2, n/2))
set.seed(0)
y = f + rnorm(n,0,sigma)

## Plot the changepoints
pt.col = 'grey50'
pt.pch = 16
lcol.mn = "red"
lwd.mn = 2
lty.cps = 2
lcol.cps = "lightgrey"
plot(y, col = pt.col, pch = pt.pch, axes=FALSE)
axis(1);axis(2)
cps = changepoints(sbs(y) , Kmax=2)$cpt.th
cps = unlist(cps)
lines(get.means(y, cps), col = lcol.mn, lwd = lwd.mn)
