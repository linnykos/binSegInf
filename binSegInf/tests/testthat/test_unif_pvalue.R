
context("Test wildBinSeg.R and some helper functions.")

test_that("Null p-values are all uniform", {

    ## Test settings

    ## Load in simulation driver functions
    source('../main/justin/sim-driver.R')

    n = 10
    numSteps = 1
    numIntervals = 100
    nsim=1000
    lev=0
    sigma=1
    nreplicate = 100
    mc.cores=3
    nsim.is=1000
    sim.settings = list(numIntervals=numIntervals,
                        nreplicate=nreplicate,
                        lev=lev,
                        numSteps=numSteps,
                        nsim.is=nsim.is)
    nsim.is = 100

    ## Simulation settings
    sim.settings = list(sigma=1, lev=0, nsim.is=10, numSteps=1,
                        numIntervals=20, n=6, meanfun=onejump,
                        reduce=FALSE,augment=TRUE,  bootstrap=FALSE, std.bootstrap=NULL,
                        cleanmn.bootstrap=NULL, thresh = 1,
                        type = "random")##plain
    sim.settings.plain = sim.settings; sim.settings.plain[["type"]]="plain"

    ## Actually run the simulations
    printprogress <- function(isim,nsim){cat("\r", "simulation ", isim,
                                             "out of", nsim)}
    methods = list(onesim_bsft, onesim_bsfs, onesim_wbs, onesim_wbs,
                   onesim_fusedlasso, onesim_fusedlasso)
    settings = list(sim.settings, sim.settings, sim.settings.plain,
                    sim.settings, sim.settings.plain, sim.settings)

    ## Conduct the ks. tests
    Map(function(mymethod, mysetting){
        a = mclapply(1:nsim,
                     function(isim){printprogress(isim,nsim); onesim_bsft(sim.settings)},
                     mc.cores=3)
        ## qqunif(unlist(a))
        expect_equal(ks.test(unlist(a),punif)$p.value<0.05, FALSE)
    }, methods, settings)


  ## Erase when done:
    a1 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_bsft(sim.settings)}, mc.cores=3)
    a2 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_bsfs(sim.settings)}, mc.cores=3)
    a3 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_wbs(sim.settings.plain)}, mc.cores=3)
    a4 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_wbs(sim.settings)}, mc.cores=3)
    a5 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_fusedlasso(sim.settings.plain)}, mc.cores=3)
    a6 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_fusedlasso(sim.settings)}, mc.cores=3)

    ## Plot and test
    methodnames = c("bsft", "bsfs", "wbs-plain", "wbs-rand", "fusedlasso-plain", "fuselasso-rand")
    for(ii in 1:6){
        ## qqunif(unlist(a))
        cat("testing", methodnames[ii])
        a = list(a1,a2,a3,a4,a5,a6)[[ii]]
        expect_equal(ks.test(unlist(a),punif)$p.value<0.05, FALSE)
    }



    ## Ideas:: Randomization wrapper to methods that produce obj$cp, obj$cp.sign? Or
    ## randomization wrapper once given a v?
