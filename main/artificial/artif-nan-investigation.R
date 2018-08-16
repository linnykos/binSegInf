# Data directory
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))


## Make a single simulation example where you will get location 1859 (or 1862;
## not sure why it would differ this time around)
set.seed(0)
sigma = sd(y.orig[1:200])
y = newmn[-(1:200)] + bootstrap_sample(resid.cleanmn[-(1:200)])

numIntervals=length(y)
consec=2
postprocess=TRUE
better.segment=TRUE
locs=1862
numIS=100; inference.type="pre-multiply"=
improve.nomass.problem=TRUE; bits=1000; write.time=TRUE
verbose=TRUE
how.close=5

    max.numSteps = 20
    g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=max.numSteps,
                              inference.type='none')
    cumsum.y = cumsum(y)

    ## Collect the IC information and polyhedron
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    ic_poly = ic_obj$poly
    stoptime  = ic_obj$stoptime

    if(verbose) cat("stoptime is", stoptime, fill=TRUE)

    ## Check for flag
    if(ic_obj$flag!="normal" ){
        return(NULL)
    }

    ## Extract changepoints from stopped model and declutter
    cp = g$cp[1:stoptime]
    cp.sign = g$cp.sign[1:stoptime]

    if(postprocess){
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
    }

    ## # ## Plot things
    ## plot(y)
    ## abline(v=cp, col='purple',lwd=2)
    ## abline(v=g$cp, col="grey80")
    ## ## abline(v=g$cp[1:ic_obj$stoptime], col='blue', lwd=2)
    ## text(x=g$cp+3, y=rep(1,length(g$cp)), label = 1:length(g$cp))


    ## Form contrasts
    if(better.segment){
        vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g, cps=cp)
    } else {
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y))
    }

    ## Retain only the changepoints we want results from:
    print((abs(as.numeric(names(vlist)))))
    retain = which((abs(as.numeric(names(vlist))) %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))

    ## Temporarily added
    if(length(retain)!=0) browser()

    ## Calculate the p-values
    vlist = vlist[retain]
iv = 1
     v = vlist[[iv]]
cumsum.v = cumsum(v)
pv = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                      sigma=sigma,
                                      numIS=numIS,
                                      inference.type=inference.type,
                                      cumsum.y=cumsum.y,
                                      cumsum.v=cumsum.v,
                                      stop.time=stoptime+consec,
                                      ic.poly=ic_poly,
                                      improve.nomass.problem=improve.nomass.problem,
                                      bits=bits))
        if(write.time) write.time.to.file(myfile="rwbs-main-example-timing.txt")
        return(pv)})
    names(pvs) = names(vlist)
