## Synopsis: Compare criteria of fused lasso and binary segmentation
source("helpers.R")
outputdir = "../figures"

## Generate date
mn = c(rep(0,5), rep(lev, 30), rep(0,5), rep(lev, 30))
set.seed(3)
y = mn + rnorm(length(mn), 0, 1)
la("~/repos/genlassoinf/genlassoinf/")
D = makeDmat(length(y),order=0)
out = dualpathSvd4(y, D=D, maxsteps=4)
print(out$cp)
crit.binseg = out$binseg.coeflist[[istep]]%*%y
print(which.max(abs(crit.binseg)/max(abs(crit.binseg))))

## Step 1
for(istep in 1:3){
    crit.binseg = out$binseg.coeflist[[istep]]%*%y
    crit.fl = out$coeflist[[istep]]%*%y
    makejpg(outputdir, paste0("step",istep,".jpg"), 500, 751)
    par(mfrow=c(2,1))
    plot(out$fusedlasso.weightlist[[istep]], type='l', lwd=2, ylab="")
    lines(out$binseg.weightlist[[istep]],col='red', type='l', lwd=2)
    crit.base = out$base.coeflist[[istep]] %*% y
    lines(abs(crit.base)/max(abs(crit.base))*10, col='blue', lwd=2)
    legend("topright", col=c("black","red","blue"), lwd=c(2,2,2),
           legend=c("Fused Lasso w_i", "Binseg w_i", "|Mean Difference|"), bg="transparent")
    title(main=paste0("Right-to-left mean difference\n and weights for each algorithm, at step ", istep))
    plot(abs(crit.binseg)/max(abs(crit.binseg))*100, type='l', col='red', lwd=2, ylab="")
    points(crit.fl/max(crit.fl)*100, type='l', col='black', lwd=2)
    title(main=paste0("Criteria at step ", istep))
    legend("topright", col=c("black","red"), lwd=c(2,2),
           legend=c("Fused Lasso criteria", "Binseg criteria"), bg="transparent")
    graphics.off()
}

## Make data plot
jpeg(file.path(outputdir, "sampledat.jpg"), width=500, height=500)
plot(y)
lines(mn, col='red', lwd=3)
graphics.off()









## fl.critfun = function(j,y){j/n*sum(y[(j+1):n]) - (n-j)/n*sum(y[1:j])}
## fl.critfun(5,y)
## fl.critfun(40,y)
## y = rep(0,n)

## nrep=10000
## crits.base = crits.binseg = crits.fl = matrix(NA,ncol=n-1, nrow=nrep)
## fac.binseg = sapply(1:(n-1),function(i){sqrt(i*(n-i)/n)})
## fac.fl = sapply(1:(n-1),function(i){i*(n-i)/n})
## for(irep in 1:nrep){
##     y = rnorm(n,0,1)
##     crits.base[irep,] = as.numeric(out.new$M%*%y)
##     crits.binseg[irep,] = as.numeric(out.new$M%*%y)*fac.binseg
##     crits.fl[irep,] = as.numeric(out.new$M%*%y)*fac.fl
## }
## cp.binseg = apply(crits.binseg,1,function(myrow){which.max(abs(myrow))})
## hist(cp.binseg)

## ## This is happening because ?
## plot(colMeans(crits.binseg), ylim = c(-1,1))
## boxplot(abs(crits.binseg))
## boxplot(abs(crits.fl))







## ## Does this match what Ryan's one-step fused lasso code is doing?

