source("simulation_header.R")

trials <- 1000
num.height <- 21
n <- 100
cores <- 20

jump.height <- exp(seq(log(0.05), log(5), length.out = num.height))
paramMat <- as.matrix(expand.grid(0, jump.height))

rule_closure <- function(n, method = binSeg_fixedSteps){
  function(vec){
    dat <- CpVector(n, vec[c(1,2,1)], c(1/3, 2/3))
    y <- dat$data
    
    obj <- method(y, 1)
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
    res <- pvalue(y, poly, contrast)
    
    obj2 <- method(y, 2)
    poly2 <- polyhedra(obj2, y)
    contrast_a <- contrast_vector(obj2, 1, sorted = T)
    contrast_b <- contrast_vector(obj2, 2, sorted = T)
    res_a <- pvalue(y, poly2, contrast_a)
    res_b <- pvalue(y, poly2, contrast_b)
    
    truth <- binSegInf:::.formMeanVec(n, dat$jump.height, dat$jump.idx)
    
    c(res, jumps(obj), sum((obj$y.fit - truth)^2)/n,
      res_a, res_b, jumps(obj2, sorted = T), sum((obj2$y.fit - truth)^2)/n)
  }
}

rule_ss_closure <- function(n, method = binSeg_fixedSteps){
  function(vec){
    dat <- CpVector(n, vec[c(1,2,1)], c(1/3, 2/3))
    y <- dat$data
    
    obj <- sample_splitting(y, method = method, numSteps = 1)
    contrast <- contrast_vector_ss(obj, 1)
    res <- pvalue_ss(y, contrast)
    
    obj2 <- sample_splitting(y, method = method, numSteps = 2)
    contrast_a <- contrast_vector_ss(obj2, 1, sorted = T)
    contrast_b <- contrast_vector_ss(obj2, 2, sorted = T)
    res_a <- pvalue_ss(y, contrast_a)
    res_b <- pvalue_ss(y, contrast_b)
    
    c(res, jumps(obj),
      res_a, res_b, jumps(obj2, sorted = T))
  }
}

############################

rule_bsFs <- rule_closure(n, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(n, method = fLasso_fixedSteps)
rule_ss <- rule_ss_closure(n)
criterion <- function(x, vec){x}

bsFs_1JumpPValue <- simulationGenerator(rule_bsFs, paramMat, criterion,
                                        trials, cores)
save.image(file = paste0("res/pvalue_oneJump_", Sys.Date(), ".RData"))

flFs_1JumpPValue <- simulationGenerator(rule_flFs, paramMat, criterion,
                                        trials, cores)
save.image(file = paste0("res/pvalue_oneJump_", Sys.Date(), ".RData"))

ss_1JumpPValue <- simulationGenerator(rule_ss, paramMat, criterion,
                                      trials, cores)
save.image(file = paste0("res/pvalue_oneJump_", Sys.Date(), ".RData"))
quit(save = "no")
