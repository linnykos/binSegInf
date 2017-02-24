source("../main/simulation_header.R")

trials <- 1000
num.height <- 21
num.loc <- 11
n <- 100
cores <- 20

jump.height <- exp(seq(log(0.05), log(5), length.out = num.height))
jump.loc <- seq(0, 1, length.out = num.loc)[2:(num.loc-1)]
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc))

rule_closure <- function(n, gridsize = 250, method = binSeg_fixedSteps){
  function(vec){
    dat <- CpVector(n, vec[1:2], vec[3])
    y <- dat$data
    
    obj <- method(y, 1)
  
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- tryCatch({
      c(confidence_interval(y, poly, contrast, gridsize = gridsize), 0)
    }, warning = function(w){
      c(-Inf, Inf, -1)
    })
    
    truth <- binSegInf:::.formMeanVec(n, dat$jump.height, dat$jump.idx)
    c(abs(contrast %*% y), abs(contrast %*% truth), res[1:2], jumps(obj),
      sum((obj$y.fit - truth)^2)/n, res[3])
  }
}

rule_ss_closure <- function(n, method = binSeg_fixedSteps){
  function(vec){
    dat <- CpVector(n, vec[1:2], vec[3])
    y <- dat$data
    
    obj <- sample_splitting(y, method = method, numSteps = 1)
    contrast <- contrast_vector_ss(obj, 1)
    
    res <- confidence_interval_ss(y, contrast)
    
    truth <- binSegInf:::.formMeanVec(n, dat$jump.height, dat$jump.idx)
    c(abs(contrast %*% y), abs(contrast %*% truth), res[1:2], jumps(obj), 
      sum((obj$y.fit)^2)/n)
  }
}

############################

rule_bsFs <- rule_closure(n, 250, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(n, 250, method = fLasso_fixedSteps)
rule_ss <- rule_ss_closure(n, method = binSeg_fixedSteps)
criterion <- function(x, vec){x}

bsFs_1JumpCI <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/CI_oneJump_", Sys.Date(), ".RData"))

flFs_1JumpCI <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/CI_oneJump_", Sys.Date(), ".RData"))

ss_1JumpCI <- simulationGenerator(rule_ss, paramMat, criterion,
                                    trials, cores)
save.image(file = paste0("../results/CI_oneJump_", Sys.Date(), ".RData"))
quit(save = "no")
