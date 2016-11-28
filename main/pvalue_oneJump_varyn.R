source("simulation_header.R")

trials <- 1000
num.n <- 21
num.loc <- 11
cores <- 20

jump.height <- 0.25
jump.loc <- seq(0, 1, length.out = num.loc)[2:(num.loc-1)]
n.vec <- round(exp(seq(log(10), log(5000), length.out = num.n)))
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc, n.vec))

rule_closure <- function(method = binSeg_fixedSteps){
  function(vec){
    y <- CpVector(vec[4], vec[1:2], vec[3])$data
    
    obj <- method(y, 1)
  
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- pvalue(y, poly, contrast)
    c(res, jumps(obj))
  }
}

############################

rule_bsFs <- rule_closure(method = binSeg_fixedSteps)
rule_flFs <- rule_closure(method = fLasso_fixedSteps)
criterion <- function(x, vec){x}

bsFs_1JumpPValue_nvary <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
flFs_1JumpPValue_nvary <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)

save.image(file = paste0("res/pvalue_oneJump_varyn_", Sys.Date(), ".RData"))
quit(save = "no")
