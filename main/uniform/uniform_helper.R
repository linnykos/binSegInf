## Helper function
make.qq.plot <- function(myresult, nsim=1000){

    ## Test code:
    nsim = 1000
    myresult = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        dosim_compare(type="sbs.nonrand", n=20, lev=0, numIS=100, meanfun=fourjump,
                      visc=1:20, numSteps=1, bits=1000,
                      max.numIS=2000)
    }, mc.cores=4)

    pvs = sapply(myresult, function(a) a$pvs)
    qqunif(pvs)

}
