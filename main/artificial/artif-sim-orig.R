## Synopsis: Run a single example on original 05296 data.
datadir = "../data"
outputdir = "../output"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))


## Estimate sigma
sigma = sd(y.orig[1:200])
y = y.orig[201:length(y.orig)]


set.seed(1)
orig.pvs.rwbs.segment = do_rwbs_inference(y=y, max.numSteps=20,
                                  numIntervals=length(y), consec=2,
                                  sigma=sigma, postprocess=FALSE,
                                  better.segment=FALSE, locs=1:length(y),
                                  numIS=100, inference.type="pre-multiply",
                                  improve.nomass.problem=TRUE, bits=1000,
                                  write.time=FALSE, verbose=TRUE)

save(orig.pvs.rwbs.segment, file=file.path(outputdir,"rwbs-orig-segment.Rdata"))


## ## Using segmentplus
## set.seed(1)
## orig.pvs.rwbs.segmentplus = do_rwbs_inference(y=y, max.numSteps=20,
##                                   numIntervals=length(y), consec=2,
##                                   sigma=sigma, postprocess=FALSE,
##                                   better.segment=TRUE, locs=1:length(y),
##                                   numIS=100, inference.type="pre-multiply",
##                                   improve.nomass.problem=TRUE, bits=1000,
##                                   write.time=FALSE, verbose=TRUE)

## save(orig.pvs.rwbs.segmentplus,
##      file=file.path(outputdir,"rwbs-orig-segmentplus.Rdata"))



## Using decluttering as well
## set.seed(1)
## orig.pvs.rwbs.declutter = do_rwbs_inference(y=y, max.numSteps=20,
##                                   numIntervals=length(y), consec=2,
##                                   sigma=sigma, postprocess=TRUE,
##                                   better.segment=TRUE, locs=1:length(y),
##                                   numIS=100, inference.type="pre-multiply",
##                                   improve.nomass.problem=TRUE, bits=1000,
##                                   write.time=FALSE, verbose=TRUE)

save(orig.pvs.rwbs.segmentplus,
     file=file.path(outputdir,"rwbs-orig-declutter.Rdata"))
