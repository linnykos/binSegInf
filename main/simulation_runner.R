rm(list=ls())

library(devtools)
install_github("linnylin92/binSegInf", ref = "kevin", subdir = "binSegInf")

library(binSegInf)
library(foreach)
library(doMC)

source("simulation_base.R")

rm(list=ls())
source("pvalue_oneJump.R")

quit(save = "no")