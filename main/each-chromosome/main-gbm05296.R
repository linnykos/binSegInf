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

## Estimate sd from these
chrome2 = chrome.data[chrome.data[,"chr"]==2,"y"]
chrome3 = chrome.data[chrome.data[,"chr"]==3,"y"]
chrome5 = chrome.data[chrome.data[,"chr"]==5,"y"]
chrome6 = chrome.data[chrome.data[,"chr"]==6,"y"]
sd(c(chrome2,chrome3, chrome5,chrome6))

## Do inference on these
chrome1 = chrome.data[chrome.data[,"chr"]==1,"y"]
chrome2 = chrome.data[chrome.data[,"chr"]==2,"y"]
chrome4 = chrome.data[chrome.data[,"chr"]==4,"y"]
chrome10 = chrome.data[chrome.data[,"chr"]==10,"y"]
chrome11 = chrome.data[chrome.data[,"chr"]==11,"y"]
four.chrome.dats = list(chrome1= chrome1, chrome4=chrome4, chrome10=chrome10, chrome11=chrome11)
header ="gm05296"
mc.cores=3

## Run main driver
source("../main/each-chromosome/gbm-driver.R")
