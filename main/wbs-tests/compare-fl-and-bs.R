## Synopsis 1: Compare FL and BS (nonrand) version with fixed four steps to
## investigate detection and conditional power.

source("../main/wbs-tests/sim-helpers.R")
## type = match.arg(type)
nsim=400
results = mclapply(1:nsim,function(isim){printprogress(isim,nsim);dosim_compare_fl_and_bs()}, mc.cores=4)

## Detection of FL seems to be better than that of BS.
detection.sbs <- sapply(results, function(myresult){    return(length(myresult$pvs.sbs.nonrand)/4) })
detection.fl <- sapply(results, function(myresult){    return(length(myresult$pvs.fl.nonrand)/4) })
print(mean(detection.sbs))
print(mean(detection.fl))

## Conditional power
cond.pow.sbs <- sapply(results, function(myresult){return(mean(myresult$pvs.sbs.nonrand < 0.05/4)) })
cond.pow.fl <- sapply(results, function(myresult){return(mean(myresult$pvs.fl.nonrand < 0.05/4))})
cond.pow.fl <- cond.pow.fl[which(!is.nan(cond.pow.fl))]
print(cond.pow.sbs)
print(cond.pow.fl)

## The conclusion is clear: with four steps, the conditional power and detection are both


## Compare with more ideal but comparable conditions for the two methods: with
## decluttering and with IC stopping:
nsim=100
results = mclapply(1:nsim, function(isim){printprogress(isim,nsim);dosim_compare_fl_and_bs_with_stoprule_and_decluttering()}, mc.cores=4)
sapply(results, function(myresult){myresult$pvs.fl.nonrand})


## Do it once and check.
source("../main/wbs-tests/sim-helpers.R")
set.seed(1)
a = dosim_compare_fl_and_bs_with_stoprule_and_decluttering()

## Why do I get NaNs? Deal with it tomorrow...

