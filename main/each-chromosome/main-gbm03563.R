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


chrome2 = chrome.data[chrome.data[,"chr"]==2,"y"] ## This is the noise-harvesting one
chrome4 = chrome.data[chrome.data[,"chr"]==4,"y"] ## This is the noise-harvesting one
chrome5 = chrome.data[chrome.data[,"chr"]==5,"y"] ## This is the noise-harvesting one
chrome6 = chrome.data[chrome.data[,"chr"]==6,"y"] ## This is the noise-harvesting one
sd(chrome2)
sd(c(chrome2, chrome4, chrome5, chrome6))


## Repeat the same analysis for four lines Try 1,3,9,11 , where 3 and 9 have
## jumps according to Olshen 2004.
chrome1 = chrome.data[chrome.data[,"chr"]==1,"y"]
chrome3 = chrome.data[chrome.data[,"chr"]==3,"y"]
chrome9 = chrome.data[chrome.data[,"chr"]==9,"y"]
chrome11 = chrome.data[chrome.data[,"chr"]==11,"y"]
header = "gm03563"
mc.cores = 4

four.chrome.dats = list(chrome1=chrome1, chrome3=chrome3,
                        chrome9=chrome9, chrome11=chrome11)

four.chrome.names = paste("Chromosome", c(1,3,9,11))

## Source in driver
source("../main/each-chromosome/gbm-driver.R")
