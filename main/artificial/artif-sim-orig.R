## Synopsis: Run a single example on original 05296 data.
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))


## Estimate sigma
sigma = sd(y.orig[1:100])
y.orig = y.orig[101:length(y.orig)]

## Pass in y
orig.pvs.rwbs = do_rwbs_inference(y=y.orig, max.numSteps=10,
                                  numIntervals=length(y.orig), consec=2,
                                  sigma=sigma, postprocess=TRUE,
                                  better.segment=TRUE, locs=1:n,
                                  numIS=100, inference.type="pre-multiply",
                                  improve.nomass.problem=TRUE, bits=1000,
                                  write.time=FALSE, verbose=TRUE)
