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

## Load results
load(file=file.path(outputdir, outputfilename))

pvtable.list = list()
## Make plot
for(ii in 1:length(four.chrome.dats)){

    ## Plot settings
    y = four.chrome.dats[[ii]]
    w=5.5; h=4
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


    ## Plot data and lines
    wbs.result = wbs.results[[ii]]
    locs = (wbs.result)$locs.retained
    pvs = (wbs.result)$pvs
    points(y, col=col.dat, pch=pch.dat)
    col.vline = rep("grey80", length(locs))
    which.sig = which(pvs<0.05/length(locs))
    col.vline[which.sig] = "grey30"
    lty.vline = rep(2, length(locs))
    lty.vline[which.sig] = 1
    abline(v=abs(locs), col=col.vline, lty=lty.vline, lwd=lwd.vline)
    lines(get_piecewise_mean(y,sort(abs(locs))), col='red', lwd=2)

    ## Add letters
    text(x=abs(locs)+1,
         y=rep(0.7, length(locs)),
         label = toupper(letters[1:length(locs)]))

    graphics.off()

    pvtable = rbind(locs, pvs)
    colnames(pvtable) = toupper(letters[1:length(locs)])
    pvtable.list[[ii]] = pvtable
}

## Make tables.
for(ii in 1:length(four.chrome.dats)){
    pvtable = pvtable.list[[ii]]
    print(names(four.chrome.dats)[ii])
    print(xtable::xtable(pvtable, col.names=FALSE))
}
save(four.chrome.dats, wbs.results, pvtable.list, file=file.path(outputdir, outputfilename))
