source("simulation_header.R")

trials <- 250
num.height <- 21
num.loc <- 11
n <- 100
cores <- 20

jump.height <- exp(seq(log(0.05), log(5), length.out = num.height))
jump.loc <- seq(0, 1, length.out = num.loc)[2:(num.loc-1)]
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc))

rule_bsFs_closure <- function(n, gridsize = 100){
  function(vec){
    dat <- CpVector(n, vec[1:2], vec[3])
    y <- dat$data
    
    obj <- binSeg_fixedSteps(y, 1)
  
    poly <- form_polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- confidence_interval(y, poly, contrast, gridsize = gridsize)
    
    truth <- binSegInf:::.formMeanVec(n, dat$jump.height, dat$jump.idx)
    c(abs(contrast %*% y), abs(contrast %*% truth), res, get_jumps(obj))
  }
}

############################

rule_bsFs <- rule_bsFs_closure(n, 50)
criterion <- function(x, vec){x}

bsFs_1JumpCI <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, NA)

save.image(file = paste0("res/CI_oneJump_bsFs_", Sys.Date(), ".RData"))
quit(save = "no")
