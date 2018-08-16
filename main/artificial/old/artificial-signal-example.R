load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
  newmn = newmn[seq(from=1,to=length(newmn),length=200)]
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
n = length(coriell_mn(1))
nsim.is = 50
numSteps = 5
numIntervals = n/2
n.levs = 2
levs=seq(from=0,to=3,length=n.levs)
nsims =rep(1,n.levs)
sim.settings <- list(levs = levs,
                     nsim.is = nsim.is,
                     numSteps = numSteps,
                     numIntervals = numIntervals,
                     n = n,
                     mn = coriell_mn,
                     nsims = nsims,
                     sigma = std,
                     seed = 0)

sim_driver(sim.settings,"artificial-example.Rdata", dir="../data")
load("../data/artificial-example.Rdata")

## Plot settings
xlab = "Location"
ylab = ""
w = 7; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-1.5,2)
mar = c(4.5,4.5,0.5,0.5)

for(pname in c("pmat.bsfs", "pmat.wbsfs")){
  for(ii in 1:2){

    ## Extract things
    cp = ((results[[ii]])[[pname]])[,"cp"]
    print(cp)
    ord.cp = order(cp)
    cp = cp[ord.cp]
    pv = (((results[[ii]])[[pname]])[,"pv"])[ord.cp]
    beta0 <- get_piecewise_mean(y, cp)
    Letters = toupper(letters[1:length(cp)])

    ## Plot setup
    ylim = c(range(y0*1.5))
    mn0 = sim.settings$mn(lev=sim.settings$levs[ii],n=sim.settings$n)
    set.seed(0)
    y0 = mn0 + rnorm(sim.settings$n, 0, sim.settings$sigma)

    ## Make plot
    filename = paste0("artif-example-",pname,"-",ii,".pdf")
    pdf(file.path("../main/figures", filename), width=w,height=h)
    par(mar=mar)
    plot(y0,ylim=ylim,xlab=xlab, col=pcol,pch=pch,axes=FALSE, ylab=ylab)
    axis(1);axis(2)
    lines(mn0, col='blue')
    lines(beta0, col='red')
    print(cp)
    abline(v=cp,col='lightgrey')
    text(x=cp, y=0.8*max(y0),label=Letters,cex=1)
    legend("bottomright", pch=c(pch,NA,NA,NA), lty=c(NA,1,1,2),lwd=c(NA,lwd,lwd,1),
           col = c(pcol,"blue","red", "black"),legend=c("Data","Mean", "Estimate","Changepoint"), bg="white")
    legend("topleft", legend=paste(Letters,signif(pv,3)), bg="white")
    graphics.off()
  }
}


## Example plot of entire thing
load("../data/coriell.Rdata")
coriell_mn <- function(lev=1,n){
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
ylim = c(range(y0*1.5))
mn0 = coriell_mn(lev=sim.settings$levs[ii],n=sim.settings$n)
set.seed(0)
y0 = mn0 + rnorm(length( coriell_mn()), 0, sim.settings$sigma)

filename = paste0("artif-example-full",".pdf")
pdf(file.path("../main/figures", filename), width=w,height=h)
par(mar=mar)
plot(y0,ylim=ylim,xlab=xlab, col=pcol,pch=pch,axes=FALSE, ylab=ylab)
axis(1);axis(2)
lines(mn0, col='blue')
legend("bottomright", pch=c(pch,NA,NA,NA), lty=c(NA,1,1,2),lwd=c(NA,lwd,lwd,1),
       col = c(pcol,"blue","red", "black"),legend=c("Data","Mean", "Estimate","Changepoint"), bg="white")
graphics.off()
