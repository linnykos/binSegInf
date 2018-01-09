## Synopsis: Analysis for GM03563, available from the CUMSEG package, in Fig 2
## of Olshen2004 CBS paper, and has truths in chromosome 3 and 9.

## Load data
library(cumSeg)
data(fibroblast)
source(file=file.path("../main/artificial/artif-helpers.R"))
source(file=file.path("../main/justin/sim-helper.R"))
outputdir = "../output"


y.orig = fibroblast$gm03563
chr = fibroblast$Chromosome

## Chrome data
chrome.data = data.frame(y=y.orig, chr=chr)
chrome.data = chrome.data[which(!is.na(chrome.data[,"y"])),]
## plot(chrome.data[,"y"])
## abline(v=c(0, cumsum(table(chrome.data[,"chr"]))))

## Repeat the same analysis for four lines
## Try 1,3,9,11 , where
chrome1 = chrome.data[chrome.data[,"chr"]==1,"y"]
chrome2 = chrome.data[chrome.data[,"chr"]==2,"y"] ## This is the noise-harvesting one
chrome3 = chrome.data[chrome.data[,"chr"]==3,"y"]
chrome9 = chrome.data[chrome.data[,"chr"]==9,"y"]
chrome11 = chrome.data[chrome.data[,"chr"]==11,"y"]
header = "gbm05363"

four.chrome.dats = list(chrome1=chrome1, chrome3=chrome3,
                        chrome9=chrome9, chrome11=chrome11)


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
    filename = paste0(header, "-", names(four.chrome.dats)[[ii]], ".pdf")
    pdf(file=file.path(outputdir, filename))
    plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab=xlab, ylab=ylab)
    axis(1); axis(2);
    points(y, col=col.dat, pch=pch.dat)

    ## Run the inference procedue
    ## sigma = sd(y.orig[1:200], na.rm=TRUE)
    sigma = sd(chrome2, na.rm=TRUE)
    set.seed(0)
    wbs.result = do_rwbs_inference(y=y, max.numSteps=10,
                                   numIntervals=2*length(y), consec=2,
                                   sigma=sigma, postprocess=TRUE,
                                   better.segment=FALSE, locs=1:length(y),
                                   numIS=100, inference.type="pre-multiply",
                                   improve.nomass.problem=TRUE, bits=1000,
                                   write.time=FALSE, verbose=TRUE)
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

