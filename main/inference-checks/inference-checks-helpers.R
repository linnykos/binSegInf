##' Make appropriate filename from simulation settings.
makesimfilename <- function(sim.settings){
    things = Map(function(a,b){paste0(a,":",b)}, names(sim.settings),
        sim.settings )

    things.pasted = paste(things, collapse="-")
    return(paste0(things.pasted, ".Rdata"))
}
