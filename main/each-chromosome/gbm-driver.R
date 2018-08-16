## Synopsis: Analyze individual chromosomes, for artificial mean for Coriell
## Cell lines 05296 and 13330. Sourced in from the settings files

## Do individual analyses
outputfilename = paste0(header, "-", "results.Rdata")
wbs.results = list()
for(ii in 1:length(four.chrome.dats)){

    ## Run the inference procedue
    y = four.chrome.dats[[ii]]
    ## sigma = sd(y.orig[1:200], na.rm=TRUE)
    if(header=="gm05296"){
        sigma = sd(c(chrome2,chrome3, chrome5,chrome6))
    } else {
        ## sigma = sd(chrome2, na.rm=TRUE)
        sigma = sd(c(chrome2, chrome4, chrome5, chrome6))
    }
    set.seed(2)
    wbs.result = do_rwbs_inference(y=y, max.numSteps=10,
                                   numIntervals = length(y), consec=2,
                                   sigma=sigma, postprocess=TRUE,
                                   better.segment=FALSE, locs=1:length(y),
                                   numIS=100, inference.type="pre-multiply",
                                   improve.nomass.problem=TRUE, bits=3000,
                                   write.time=FALSE, verbose=TRUE,
                                   max.numIS = 20000,
                                   mc.cores=mc.cores)
    wbs.results[[ii]] = wbs.result
    save(four.chrome.dats, wbs.results, file=file.path(outputdir, outputfilename))
}
