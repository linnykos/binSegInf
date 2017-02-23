source("simulation_header.R")

trials <- 1000
n <- 100
cores <- NA

paramMat <- as.matrix(0)

rule_closure <- function(n, gridsize = 100, method = binSeg_fixedSteps){
  function(void){
    dat <- CpVector(n, 0, NA)
    y <- dat$data
   
    obj <- method(y, 1)
  
    poly <- polyhedra(obj)
    contrast <- contrast_vector(obj, 1)

    res <- tryCatch({
      c(confidence_interval(y, poly, contrast, gridsize = gridsize), 0)
    }, warning = function(w){
      c(-Inf, Inf, -1)
    })
    
    c(abs(contrast %*% y), res[1:2], jumps(obj), sum((obj$y.fit)^2)/n, res[3])
  }
}

rule_ss_closure <- function(n, method = binSeg_fixedSteps){
  function(void){
    dat <- CpVector(n, 0, NA)
    y <- dat$data
    
    obj <- sample_splitting(y, method = method, numSteps = 1)
    contrast <- contrast_vector_ss(obj, 1)
    
    res <- confidence_interval_ss(y, contrast)
    
    c(abs(contrast %*% y), res[1:2], jumps(obj))
  }
}

############################

rule_bsFs <- rule_closure(n, 50, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(n, 50, method = fLasso_fixedSteps)
rule_ss <- rule_ss_closure(n, method = binSeg_fixedSteps)
criterion <- function(x, vec){x}

bsFs_0JumpCI <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/CI_noJump_", Sys.Date(), ".RData"))

flFs_0JumpCI <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)
save.image(file = paste0("../results/CI_noJump_", Sys.Date(), ".RData"))

ss_0JumpCI <- simulationGenerator(rule_ss, paramMat, criterion,
                                    trials, cores)
save.image(file = paste0("../results/CI_noJump_", Sys.Date(), ".RData"))
quit(save = "no")
