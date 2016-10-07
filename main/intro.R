library(wbs)
library(binSegInf)

## Generate data
sigma = 1
lev1 = 0
lev2 = 3
n=100
f = c(rep(lev1,n/2), rep(lev2, n/2))
set.seed(5)
y = f + rnorm(n,0,sigma)

## Plot the changepoints
pt.col = 'grey50'
pt.pch = 16
lcol.mn = "red"
lwd.mn = 2
lty.cps = 2
lcol.cps = "lightgrey"
plot(y, col = pt.col, pch = pt.pch, axes=F)
axis(1);axis(2)
cps = changepoints(sbs(y) , Kmax=2)$cpt.th
cps = unlist(cps)
lines(get.means(y, cps), col = lcol.mn, lwd = lwd.mn)
abline(v = cps, lty = lty.cps, col = lcol.cps)


## Get the permutation t-test p-values
spurious.cp = sort(cps)[1]
real.cp = sort(cps)[2]
vec1 = y[1:spurious.cp]
vec2 = y[(spurious.cp+1):real.cp]
signif(perm.t.test(vec1,vec2,1000),3)

## Get the selective test p-values.



## Plot them.
outputdir = "figures"
filename = "abc"
pdf(filename)
