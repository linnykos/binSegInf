## Synopsis: Script to run the entire set of experiments. Run by ``Rscript
## ../compare/compare-run.R 1 2 3 4 5 6 7 8 9 10 11 12'', from the directory
## binSegInf/binSegInf

## Load the current package, helpers, and data.
load_all() 
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))
source(file=file.path("../main/bootstrap/bootstrap.R"))
outputdir="../output"


## Simulation settings
bits = 5000
mc.cores = 6
args = commandArgs(trailingOnly=TRUE)
facs = as.numeric(args)
nsims = seq(from=2000,to=4000, length=5)
