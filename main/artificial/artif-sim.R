load("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1, n){
    ## newmn = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}

## simulation settings
n = length(coriell_mn(lev=1))
numIS = 100
numSteps = 1 ## 5
numIntervals = n/2
n.levs = 1

## Generate nosie around the data
sigma=.1
y = newmn + rnorm(n, 0, sigma)
cumsum.y = cumsum(y)

## Fit initial WBS for a generous number of steps
numSteps=10
numIntervals=n
g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps,
                          inference.type='none')

## Collect the IC information and polyhedron
## Get ic-stopping polyhedron
consec=2
ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
## ic_poly = ic_to_poly(ic_obj)
ic_poly = ic_obj$poly

## Check for flag
if(ic_obj$flag=="normal" ){

    if(!randomized){
        ## Get ic-stopped model selection polyhedron
        stopped.gamma = do.call(rbind, g$rows.list[1:(ic_obj$stoptime+consec)])
        stopped.u = rep(0, nrow(stopped.gamma))
        poly = polyhedra(obj=rbind(stopped.gamma, ic_obj$poly$gamma),
                         u=c(stopped.u, ic_obj$gamma$u))
    }
} else {
    return(ic_obj$flag)
}
stoptime  = ic_obj$stoptime

## Extract changepoints from stopped model and form contrasts
cp = g$cp[1:stoptime]
cp.sign = g$cp.sign[1:stoptime]
if(better.segment){
    vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g)
} else {
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
}

## Retain only the changepoints we want results from:
retain = which((cp %in% locs))
if(length(retain)==0) return(list(pvs=c(), null.true=c()))

## Calculate the p-values
vlist = vlist[retain]
pvs = sapply(vlist, function(v){
    if(randomized){
        cumsum.v = cumsum(v)
        return(suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                                sigma=sigma,
                                                numIS=numIS,
                                                inference.type=inference.type,
                                                cumsum.y=cumsum.y,
                                                cumsum.v=cumsum.v,
                                                stop.time=stoptime+consec,
                                                ic.poly=ic_poly,
                                                improve.nomass.problem=improve.nomass.problem)
                                                ))





levs=seq(from=0,to=3,length=n.levs)
nsims = seq(from=100,to=50,length=n.levs)
bootstrap=TRUE
reduce=TRUE
mc.scores = 1




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
