load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    ## newmn = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}

n = length(coriell_mn(1))
nsim.is = 1#500
numSteps = 5
numIntervals = 500
n.levs = 1
levs=seq(from=0,to=3,length=n.levs)
nsims = seq(from=100,to=50,length=n.levs)
bootstrap=TRUE
reduce=TRUE
mc.scores = 1

sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = coriell_mn,
                     nsims = nsims,
                     sigma = std,
                     std = std,
                     bootstrap=bootstrap,
                     resid.cleanmn = resid.cleanmn
                     )


sim_driver(sim.settings = sim.settings,
           ## filenames = paste0("artificial-lev-",mylev,".Rdata"),
           filenames = "artificial.Rdata",
           dir="../data", reduce=reduce )


## Plot things
pdf("~/Desktop/sample-data.pdf",width=10,height=10)
par(mfrow=c(2,2))
for(lev in 1:4){
mn <- coriell_mn
set.seed(0)
y <- mn(lev,n) + rnorm(n,0,std)
plot(y,ylim=c(-1,1), main = paste0("Signal size /stretched/ from snr=1 to ",lev),pch=16,cex=.5,col='grey50')
lines(mn(lev,n),col='red')
}
graphics.off()
