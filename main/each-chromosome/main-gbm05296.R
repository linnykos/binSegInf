## Synopsis: Analyze individual chromosomes, for artificial mean for Coriell
## Cell lines 05296 and 13330.
library(DNAcopy)
library(genlasso); library(genlassoinf); ##library(tidyverse)
data(coriell)
datadir = "../data"
outputdir = "../output"
source(file=file.path("../main/artificial/artif-helpers.R"))
source(file=file.path("../main/justin/sim-helper.R"))

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23
CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))
y.orig = (coriell[,"Coriell.05296"])

chrome.orig = coriell[,"Chromosome"]

## Get rid of outlier
outlier.ind = which(y.orig[1:1000] < -0.5)
y.orig = y.orig[-outlier.ind]
chrome.orig = chrome.orig[-outlier.ind]

## Chrome data
chrome.data = data.frame(y=y.orig, chr=chrome.orig)
chrome.data = chrome.data[which(!is.na(chrome.data[,"y"])),]
## plot(chrome.data[,"y"])
## abline(v=c(0, cumsum(table(chrome.data[,"chr"]))))

chrome1 = chrome.data[chrome.data[,"chr"]==1,"y"]
chrome2 = chrome.data[chrome.data[,"chr"]==2,"y"]
chrome4 = chrome.data[chrome.data[,"chr"]==4,"y"]
chrome10 = chrome.data[chrome.data[,"chr"]==10,"y"]
chrome11 = chrome.data[chrome.data[,"chr"]==11,"y"]
four.chrome.dats = list(chrome1= chrome1, chrome4=chrome4, chrome10=chrome10, chrome11=chrome11)

## Do individual analyses
wbs.results = list()
for(ii in 1:length(four.chrome.dats)){
    print(names(four.chrome.dats)[ii])
    y = four.chrome.dats[[ii]]

    ## Plot settings
    w=7; h=5
    xlim = c(0, length(y))
    ylim = range(y)
    col.dat = "grey50"
    pch.dat = 16
    xlab = "Coordinate"
    ylab = expression(log ~ 2 ~ ratio)
    col.vline="grey80"
    lty.vline=2
    lwd.vline=2

    ## Make basic plot
    filename = paste0("gm05296-", names(four.chrome.dats)[[ii]], ".pdf")
    pdf(file=file.path(outputdir, filename))
    plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab=xlab, ylab=ylab)
    axis(1); axis(2);
    points(y, col=col.dat, pch=pch.dat)

    ## Run the inference procedue
    ## sigma = sd(y.orig[1:200], na.rm=TRUE)
    sigma = sd(chrome2, na.rm=TRUE)
    set.seed(2)
    wbs.result = do_rwbs_inference(y=y, max.numSteps=10,
                                   numIntervals = length(y), consec=2,
                                   sigma=sigma, postprocess=TRUE,
                                   better.segment=FALSE, locs=1:length(y),
                                   numIS=100, inference.type="pre-multiply",
                                   improve.nomass.problem=TRUE, bits=3000,
                                   write.time=FALSE, verbose=TRUE,
                                   max.numIS = 20000)
    wbs.results[[ii]] = wbs.results

    ## Plot it
    locs = (wbs.result)$locs.retained
    pvs = (wbs.result)$pvs
    abline(v=abs(locs), col=col.vline, lty=lty.vline, lwd=lwd.vline)

    ## Plot piecewise means
    lines(get_piecewise_mean(y,sort(abs(locs))), col='red', lwd=2)

    ## Add text
    text(x=abs(locs)+1,
         y=rep(0.2, length(locs)),
         label = signif(pvs,3))
    graphics.off()
}
