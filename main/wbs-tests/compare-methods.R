## Synopsis: See randomized tests power compared to the nonrandomized version
library(genlassoinf)

## Generate p-values
dosim <- function(type=c("wbs","fl","sbs"), n, lev, numIntervals=n, sigma.add=0.2){

    type = match.arg(type)
    numSteps=1
    sigma=1
    mn = c(rep(0,n/2), rep(lev,n/2))
    y = mn + rnorm(n, 0, sigma)
        numIS = 30

    if(type=="wbs"){
        ## Fit WBS, test first jump
        g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
        poly.wbs = polyhedra(obj=g$gamma, u=g$u)
        vlist <- make_all_segment_contrasts(g)
        v.wbs = vlist[[1]]
        return(data.frame(
            pv.wbs.rand =  suppressWarnings(randomize_wbsfs(v=v.wbs, winning.wbs.obj=g, sigma=sigma, numIS=numIS)),
            pv.wbs.nonrand = poly.pval2(y=y, poly=poly.wbs, v=v.wbs, sigma=sigma)$pv,
            loc.wbs = g$cp * g$cp.sign))
    }

    if(type=="fl"){
        ## Fit FL, test first jump
        D = genlassoinf::makeDmat(n,type='tf',ord=0)
        f = genlassoinf::dualpathSvd2(y, D=D, maxsteps=1, approx=T)
        Gobj.naive = genlassoinf::getGammat.naive(obj=f, y=y, condition.step=1)
        poly.fl = polyhedra(obj=Gobj.naive$G, u=Gobj.naive$u)
        vlist <- make_all_segment_contrasts(f)
        v.fl = vlist[[1]]

            return(data.frame(pv.fl.rand = randomize_genlasso(v=v.fl, pathobj=f, sigma=sigma,
                                                              numIS=numIS, sigma.add=sigma.add, orig.poly=poly.fl),
                              pv.fl.nonrand = poly.pval2(y=y, poly=poly.fl, v=v.fl, sigma=sigma)$pv,
                              loc.fl = f$cp * f$cp.sign))
    }

    if(type=="sbs"){
        ## Fit SBS, test first jump
        h = binSeg_fixedSteps(y, numSteps=numSteps)
        poly.bs = polyhedra(h)
        vlist <- make_all_segment_contrasts(h)
        v.sbs = vlist[[1]]

            return(data.frame(pv.sbs.rand = randomize_genlasso(v=v.sbs, pathobj=h,
                                                               sigma=sigma,
                                                               numIS=numIS,
                                                               sigma.add=sigma.add,
                                                               orig.poly=poly.bs),
                              pv.sbs.nonrand = poly.pval2(y=y, poly=poly.bs, v=v.sbs, sigma=sigma)$pv,
                              loc.sbs = h$cp * h$cp.sign))
    }
}

nsim=100
results.wbs = mclapply(1:nsim, function(isim){printprogress(isim,nsim);dosim(type="wbs", n=20, lev=0)}, mc.cores=4)
nsim=300
results.fl =  mclapply(1:nsim, function(isim) {printprogress(isim,nsim);dosim(type="fl", n=60, lev=0)}, mc.cores=4)
nsim=100
results.sbs = mclapply(1:nsim, function(isim){printprogress(isim,nsim);dosim(type="sbs", n=40, lev=0)}, mc.cores=4)
qqunif(results.sbs[,2])

a = do.call(rbind, results.sbs)
a = do.call(rbind, results.wbs)
a = do.call(rbind, results.fl)
qqunif(a[which(a[,"loc.wbs"]==30),"pv.wbs.rand"])
qqunif(a[,"pv.wbs.rand"])


## ## I want locations too!
## nsim=500
## lev0 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=0, mc.cores=7)
## lev1 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=1, mc.cores=7)
## lev2 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=2, mc.cores=7)
## lev3 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7)
## ## lev5 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=5, mc.cores=7)


## save(list=c("lev0","lev1"), file="lev-part1.Rdata")
## save(list=c("lev2","lev3"), file="lev-part2.Rdata")



load("~/Desktop/lev-part1.Rdata")
load("~/Desktop/lev-part2.Rdata")

## Investigating speed
Rprof("a.out")
## pv.fl.rand = randomize_genlasso(v=v, pathobj=f, sigma=sigma, numIS=numIS, sigma.add=0.2, orig.poly=poly.fl)
pv.wbs.rand =  suppressWarnings(randomize_wbsfs(v=v.wbs, winning.wbs.obj=g, sigma=sigma, numIS=numIS))
Rprof(NULL)
summaryRprof("a.out")


## Cleaning helper function
myclean <- function(mylev, returnlocs=FALSE){
   mylev = do.call(rbind,mylev) ## Still getting some NaNs.. not sure why!
   which.loc.cols = which(names(mylev)%in%c("loc.wbs", "loc.fl", "loc.sbs"))
   mylocs = mylev[,which.loc.cols]
   mylev = mylev[,-which.loc.cols]
   mynames = colnames(mylev)
   mylev = lapply(1:ncol(mylev),function(icol)mylev[,icol])
   names(mylev) = mynames
   if(returnlocs) return(mylocs)
   return(mylev)
}


## Visualize 1: compare all the results
w=h=5
for(ilev in 1:4){
    mylev = list(lev0,lev1,lev2,lev3)[[ilev]]
    pdf(file.path("~/Desktop/", paste0(ilev,".pdf")), width=w,height=h)
    myclean(mylev)
    qqunif(myclean(mylev),cols=1:6)
    graphics.off()
}

range(myclean(lev1)[[4]])
range(myclean(lev2)[[4]])
range(myclean(lev3)[[4]])
range(myclean(lev3)[[4]])

## Visualize 2: Phase transition of power loss with more and more randomization?
## /Doesn't/ occur for WBS since /more intervals just get closer to
## comprehensive/? Confirm this.


## Visualize: See the effect of these56gg
lev0 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=0, numIntervals=2, mc.cores=7)
lev0 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=0, numIntervals=4, mc.cores=7)


## Investigate power: See the power transition from small to large additive noise
nsim=500
ns0 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=0.2)
ns1 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=0.5)
ns2 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=1)
ns3 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=1.5)
ns4 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=3)
ns5 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=5)
ns6 = mclapply(1:nsim, comparesim, nsim=nsim, n=8, lev=3, mc.cores=7, sigma.add=10)
save(list=c("ns0", "ns1", "ns2", "ns3","ns4","ns5","ns6"), file="ns.Rdata")
load(file="~/Desktop/ns.Rdata")

## Visualize it (separately)
w = h = 5
for(i.ns in 1:6){
    my.ns = list(ns0, ns1, ns2, ns3, ns4, ns5, ns6)[[i.ns]]
    pdf(file.path("~/Desktop/", paste0("ns-",i.ns,".pdf")), width=w,height=h)
    ttl = paste("additive sd =", c(0.2,0.5,1,1.5,3,5,10))[i.ns]
    qqunif(myclean(my.ns),cols=1:2, main=ttl)
    mypv = myclean(my.ns)[["pv.sbs.rand"]]
    print(sum(mypv<0.05))
    graphics.off()
}

## Visualize it.
w = h = 5
pdf(file.path("~/Desktop/", paste0("ns-rand.pdf")), width=w,height=h)
for(i.ns in 1:6){
    my.ns = list(ns0, ns1, ns2, ns3, ns4, ns5, ns6)[[i.ns]]
    mypv = myclean(my.ns)[[1]]
    pts = qqunif(mypv, plot.it=FALSE)
    if(i.ns==1){
        plot(pts,col=i.ns, cex=.2)
    } else {
        points(pts,col=i.ns, cex=.2)
    }
}
legend("bottomright", col=1:6, legend = paste("additive sd =", c(0.2,0.5,1,1.5,3.5,10)), pch = rep(1,6))
graphics.off()
