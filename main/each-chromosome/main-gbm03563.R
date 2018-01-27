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

## Source in driver
source("../main/each-chromosome/gbm-driver.R")
