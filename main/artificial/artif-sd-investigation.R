## Synopsis: Is sd stable? (short answer: after eliminating outliers, yes)
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))

## SD on raw data.
sigmas = c()
starts = seq(from=0,to=500, by=20)
for(start in starts){
    sigmas = c(sigmas, sd(y.orig[start+(1:200)]))
}
par(mfrow=c(2,1))
plot(y.orig)
lines(sigmas~starts, type='o', col='red',lwd=2)

## Eliminate outliers and see that sd is much stabler
outlier.ind = which(y.orig[1:700] < -0.5)
y.orig = y.orig[-outlier.ind]

sigmas = c()
starts = seq(from=0,to=500, by=20)
for(start in starts){
    sigmas = c(sigmas, sd(y.orig[start+(1:200)]))
}
plot(y.orig)
lines(sigmas~starts, type='o', col='red',lwd=2)
