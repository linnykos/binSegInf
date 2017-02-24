source("../main/simulation_header.R")

trials <- 150
num.n <- 21
num.loc <- 11
cores <- 20

jump.height <- 0.25
jump.loc <- seq(0, 1, length.out = num.loc)[2:(num.loc-1)]
n.vec <- round(exp(seq(log(10), log(1000), length.out = num.n)))
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc, n.vec))

rule_closure <- function(gridsize = 250, method = binSeg_fixedSteps){
  function(vec){
    n <- vec[4]
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

rule_ss_closure <- function(method = binSeg_fixedSteps){
  function(vec){
    n <- vec[4]
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

rule_bsFs <- rule_closure(250, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(250, method = fLasso_fixedSteps)
rule_ss <- rule_ss_closure(method = binSeg_fixedSteps)
criterion <- function(x, vec){x}

bsFs_1JumpCI_nvary <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/CI_oneJump_vary_n_", Sys.Date(), ".RData"))

flFs_1JumpCI_nvary <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/CI_oneJump_varyn_", Sys.Date(), ".RData"))

ss_1JumpCI_nvary <- simulationGenerator(rule_ss, paramMat, criterion,
                                          trials, cores)
save.image(file = paste0("../results/CI_oneJump_varyn_", Sys.Date(), ".RData"))
quit(save = "no")
