##  recovery comparison
onesim_recovery <- function(lev){

    ## Generate data
    n=60
    meanfun=fourjump
    mn = meanfun(lev,n)
    sigma=1
    y = mn + rnorm(n, 0, sigma)
    sigma.add = 0.2
    numSteps = 4
    new.noise = rnorm(n,0,sigma.add)
    locs = sapply(1:4, function(ii){ii/5*n + c(-1,0,1)})

    ## Get four-step sbs model
    g.bs = binSeg_fixedSteps(y + new.noise, numSteps=numSteps)
    correct1 = sum(g.bs$cp%in%locs)/length(g.bs$cp)
    error1 = (length(g.bs$cp) - sum(g.bs$cp%in%locs))/length(g.bs$cp)

    ## Get stopped wbs model
    g.wbs = wildBinSeg_fixedSteps(y + new.noise, numSteps=numSteps, numIntervals=n)
    correct2 = sum(g.wbs$cp%in%locs)/length(g.wbs$cp)
    error2 = (length(g.wbs$cp) - sum(g.wbs$cp%in%locs))/length(g.wbs$cp)

    ## Get stopped cbs model
    g.cbs = circularBinSeg_fixedSteps(y + new.noise, numSteps=numSteps/2)
    correct3 = sum(g.cbs$cp%in%locs)/length(g.cbs$cp)
    error3 = (length(g.cbs$cp) - sum(g.cbs$cp%in%locs))/length(g.cbs$cp)

    return(data.frame(correct1=correct1,error1=error1,
                      correct2=correct2,error2=error2,
                      correct3=correct3,error3=error3))
}

