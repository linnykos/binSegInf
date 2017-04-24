yload("../data/coriell.Rdata")

##' Takes a given mean and multiply the maximum to have noise*lev maximum
##' height.
coriell_mn <- function(lev=1,n){
    ## newmn2 = (newmn[1101:1300][seq(from=1,to=200,length=100)])
    h = max(abs(newmn))
    return((newmn / h * std) * lev)
}
n = length(coriell_mn(1))
nsim.is = 100
numSteps = 1
numIntervals = 5000
lev=3
nreplicate =  c(100,50,10)
nsims.is = c(10,100,1000)

ptm <- proc.time()


for(ii in 2){
  nsim.is = nsims.is[ii]
  cat("nsim.is is", nsim.is, "out of", nsims.is, fill=TRUE)
  
  set.seed(0)
  y <- coriell_mn(lev,n) + rnorm(n,0,std)
  p.wsbs.list = list()

  for(irep in 1:(nreplicate[ii])){

    cat("replicate", irep, "out of", nreplicate[ii], fill=TRUE)
    print("this replicate starts")
    print(proc.time() - ptm)

    method <- wildBinSeg_fixedSteps
    intervals <- generate_intervals(n=length(y),numIntervals=numIntervals, seed=0)
    obj <- method(y, numSteps=numSteps, intervals=intervals)

    print("ran method")
    print(proc.time() - ptm)

    poly <- polyhedra(obj)

    print("collected polyhedron")
    print(proc.time() - ptm)

    contrast <- make_all_segment_contrasts(obj)
    for(ii in 1:length(obj$cp)){
        p.wbsfs[ii] <- randomized_wildBinSeg_pv(y=y,
                                                v=contrast[[ii]], sigma=sigma,
                                                numSteps=numSteps,
                                                numIntervals=numIntervals,
                                                nsim.is=nsim.is, bits=100)
    }
    p.wbsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.wbsfs)
    colnames(p.wbsfs.nonrand) = colnames(p.wbsfs) = c("isim","cp","pv")
    p.wbsfs.list[[irep]] = p.wbsfs
  }

  ## Save each result
  filename=paste0("artificial-nsimis-",nsim.is ,".Rdata")
  save(p.wsfs.list, filename, dir="../main/data")
}
