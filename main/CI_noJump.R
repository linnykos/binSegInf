source("simulation_header.R")

trials <- 1000
n <- 100
cores <- 20

paramMat <- as.matrix(0)

rule_closure <- function(n, gridsize = 100, method = binSeg_fixedSteps){
  function(void){
    dat <- CpVector(n, 0, NA)
    y <- dat$data
    
    obj <- method(y, 1)
  
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- tryCatch({
      c(confidence_interval(y, poly, contrast, gridsize = gridsize), 0)
    }, warning = function(w){
      c(-Inf, Inf, -1)
    })
    
    c(abs(contrast %*% y), res[1:2], jumps(obj), sum((obj$y.fit)^2)/n, res[3])
  }
}

############################

rule_bsFs <- rule_closure(n, 50, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(n, 50, method = fLasso_fixedSteps)
criterion <- function(x, vec){x}

bsFs_0JumpCI <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
flFs_0JumpCI <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)


save.image(file = paste0("res/CI_noJump_", Sys.Date(), ".RData"))
quit(save = "no")
