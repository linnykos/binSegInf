source("simulation_header.R")

trials <- 1000
num.height <- 21
num.loc <- 11
cores <- 20

jump.height <- 0.25
jump.loc <- 0.5
n.vec <- exp(seq(log(10), log(5000), length.out = 20))
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc, round(n.vec)))

rule_bsFs_closure <- function(){
  function(vec){
    y <- CpVector(vec[4], vec[1:2], vec[3])$data
    
    obj <- binSeg_fixedSteps(y, 1)
  
    poly <- form_polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- pvalue(y, poly, contrast)
    c(res, get_jumps(obj))
  }
}

############################

rule_bsFs <- rule_bsFs_closure()
criterion <- function(x, vec){x}

bsFs_1JumpPValue_nvary <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)

save.image(file = paste0("res/pvalue_oneJump_varyn_bsFs_", Sys.Date(), ".RData"))
quit(save = "no")
