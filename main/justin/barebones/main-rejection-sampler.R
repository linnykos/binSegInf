## Draw data and fit model -- all /original/ information
nsim = 1500
numIntervals=4
pvs4 = c()
for(isim in 1:nsim){

    cat("simulation", isim, "out of", nsim, fill=TRUE)
    n=4
    sigma=1
    ## set.seed(1)
    lev=0
    y = c(rep(0,n/2), rep(lev,n/2)) + rnorm(n,0,sigma)
    obj <- wildBinSeg_fixedSteps(y, numSteps=1, numIntervals=numIntervals)
    cp = obj$cp * obj$cp.sign
    v <- make_all_segment_contrasts(obj)[[1]]
    v = v/sqrt(sum(v*v))


    orth.proj = diag(1,n) - v%*%t(v)
    orth.comp = orth.proj%*%y
    nrep = 300
    enough.samples = FALSE
    rej.samp.vz = c()

    while(!enough.samples){
        new.rej.samp.vz = unlist(mclapply(1:nrep, function(irep){

            ## Retain if new.z is larger than observed v'y
            printprogress(irep, nrep)
            Pvz = rnorm(n=n, mean=0, sd=sigma)
            new.z = Pvz + orth.comp
            new.obj <- wildBinSeg_fixedSteps(new.z, numSteps=1, numIntervals=numIntervals)
            if(new.obj$cp * new.obj$cp.sign == cp){
                return(v %*% new.z)
            } else {
                return(NA)
            }
        }, mc.cores=6))
        new.rej.samp.vz = new.rej.samp.vz[which(!is.na(new.rej.samp.vz))]
        rej.samp.vz = c(rej.samp.vz, new.rej.samp.vz)
        enough.samples = (length(rej.samp.vz) > 100)
    }
    pv = sum(rej.samp.vz > as.numeric(v%*%y))/length(rej.samp.vz)
    pvs4[isim] = pv
    cat(fill=TRUE)
}
save(list="pvs4", file="../output/reject-binseg.Rdata")
