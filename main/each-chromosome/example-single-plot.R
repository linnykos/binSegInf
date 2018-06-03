## Synopsis: Plot an individual chromosomes, for artificial mean for Coriell
## Cell lines 05296 and 13330. Sourced in from the settings files.

## Load results
outputfilename = paste0(header, "-", "results.Rdata")
load(file=file.path(outputdir, outputfilename))

## pvtable.list = list()
## ## Make plot
## for(ii in 1:length(four.chrome.dats)){
    ii=1

    ## Plot settings
    y = four.chrome.dats[[ii]]
    w=5.5; h=3.5
    xlim = c(0, length(y))
    ylim = range(unlist(four.chrome.dats)) + c(0,0.2)##range(y)
    col.dat = "grey50"
    pch.dat = 16
    xlab = "Coordinate"
    ylab = expression(log ~ 2 ~ ratio)
    col.vline="grey80"
    lty.vline=2
    lwd.vline=2

    ## Set up plot
    plotfilename = paste0(header, "-", names(four.chrome.dats)[[ii]], ".pdf")
    cat("Made plot in:", file.path(outputdir,plotfilename), fill=TRUE)
    pdf(file=file.path(outputdir, plotfilename), width=w, height=h)
    par(mar=c(3,2,2,1))
    plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab=xlab, ylab=ylab)
    axis(1); axis(2);
