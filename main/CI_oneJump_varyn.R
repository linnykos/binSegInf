source("simulation_header.R")

trials <- 1000
num.n <- 21
num.loc <- 11
cores <- 20

jump.height <- 0.25
jump.loc <- seq(0, 1, length.out = num.loc)[2:(num.loc-1)]
n.vec <- round(exp(seq(log(10), log(5000), length.out = num.n)))
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc, n.vec))

rule_closure <- function(gridsize = 100, method = binSeg_fixedSteps){
  function(vec){
    dat <- CpVector(vec[4], vec[1:2], vec[3])
    y <- dat$data
    
    obj <- method(y, 1)
  
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- tryCatch({
      c(confidence_interval(y, poly, contrast, gridsize = gridsize), 0)
    }, warning = function(w){
      c(-Inf, Inf, -1)
    })
    
    truth <- binSegInf:::.formMeanVec(vec[4], dat$jump.height, dat$jump.idx)
    c(abs(contrast %*% y), abs(contrast %*% truth), res[1:2], jumps(obj),
      res[3])
  }
}

############################

rule_bsFs <- rule_closure(50, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(50, method = fLasso_fixedSteps)
criterion <- function(x, vec){x}

bsFs_1JumpCI <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
flFs_1JumpCI <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)

save.image(file = paste0("res/CI_oneJump_vary_n_", Sys.Date(), ".RData"))
quit(save = "no")
