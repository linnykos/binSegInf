## Clean data (from sheet 3 of ng754-s11 from doi:10.1038/ng754
rawdat = read.csv("../data/ng754-s11-sheet3.csv")
chromnums = unique(rawdat[,"Chromosome"])

## Store separate chromosomes
alldat = lapply(chromnums,
                function(ii) rawdat[rawdat[,"Chromosome"]==ii,])
names(alldat) = paste("Chromosome", chromnums)

## Store the entire raw thing as the last entry
alldat[[length(alldat)+1]] = rawdat
names(alldat)[[length(alldat)]] = "All Chromosomes"

## Store it
save(list=c("alldat"), file="../data/ng754-s11-sheet3.Rdata")



## Example plot
load(file="./data/ng754-s11-sheet3.Rdata")

## Plot settings
w = h = 5
lcol.fl = "red"
lty.fl = 1
lwd.fl = 2
xlab = "Genome Position"
ylab = "log 2 ratio"
pcol.dat = 'grey50'
pch.dat = 16

## Make basic data plots.
for(chromenum in c(4,10)){
    chromename = paste0("Chromosome ", chromenum)

    ## Load data
    dat = alldat[[chromename]]
    x = dat[,"Position.Genome..kb"]
    y = dat[,"Log2Ratio"]
    x = x[!is.na(y)]
    y = y[!is.na(y)]

    filename = paste0("ng754-", chromenum, "-justdata.pdf")
    pdf(file.path("~/Desktop",filename),width=w,height=h)
    par(mar=c(3,2,2,2))
    plot(y ~ x, pch = pch.dat, col = pcol.dat, axes = FALSE, ylab = ylab,
         xlab = xlab, ylim=c(-1,1))
    axis(1);axis(2)
    graphics.off()
}
