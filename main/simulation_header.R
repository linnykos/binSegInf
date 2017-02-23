rm(list=ls())

library(devtools)
install_github("linnylin92/binSegInf", ref = "kevin", subdir = "binSegInf",
  force = T)

library(binSegInf); library(foreach); library(doMC)
source("../main/simulation_base.R")
