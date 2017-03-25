source("../main/simulation_header.R")

trials <- 1000
num.n <- 21
num.loc <- 11
cores <- 20

jump.height <- 0.25
jump.loc <- seq(0, 1, length.out = num.loc)[2:(num.loc-1)]
n.vec <- round(exp(seq(log(10), log(1000), length.out = num.n)))
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc, n.vec))

rule_closure <- function(method = binSeg_fixedSteps){
  function(vec){
    n <- vec[4]
    dat <- CpVector(n, vec[1:2], vec[3])
    y <- dat$data
    
    obj <- method(y, 1)
  
    poly <- polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
  
    res <- pvalue(y, poly, contrast)
    
    truth <- binSegInf:::.formMeanVec(n, dat$jump.height, dat$jump.idx)
    c(res, jumps(obj), sum((obj$y.fit - truth)^2)/n)
  }
}

rule_ss_closure <- function(method = binSeg_fixedSteps){
  function(vec){
    n <- vec[4]
    dat <- CpVector(n, vec[1:2], vec[3])
    y <- dat$data
    
    obj <- sample_splitting(y, method = method, numSteps = 1)
    contrast <- contrast_vector_ss(obj, 1)
    
    res <- pvalue_ss(y, contrast)
    
    c(res, jumps(obj), sum((obj$y.fit)^2)/n)
  }
}

############################

rule_bsFs <- rule_closure(method = binSeg_fixedSteps)
rule_flFs <- rule_closure(method = fLasso_fixedSteps)
rule_ss <- rule_ss_closure(method = binSeg_fixedSteps)
criterion <- function(x, vec){x}

bsFs_1JumpPValue_nvary <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/pvalue_oneJump_varyn_", Sys.Date(), ".RData"))

flFs_1JumpPValue_nvary <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/pvalue_oneJump_varyn_", Sys.Date(), ".RData"))

ss_1JumpPValue_nvary <- simulationGenerator(rule_ss, paramMat, criterion,
                                      trials, cores)
save.image(file = paste0("../results/pvalue_oneJump_varyn_", Sys.Date(), ".RData"))


quit(save = "no")
