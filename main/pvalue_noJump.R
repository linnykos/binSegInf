source("simulation_header.R")

trials <- 1000
n <- 100
cores <- 20

paramMat <- as.matrix(0)

rule_closure <- function(n, method = binSeg_fixedSteps){
  function(void){
    dat <- CpVector(n, 0, NA)
    y <- dat$data
    
    obj <- method(y, 1)
  
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
    
    if(any(poly$gamma %*% y < poly$u)) stop()
  
    res <- pvalue(y, poly, contrast)
    c(res, jumps(obj), sum((obj$y.fit)^2)/n)
  }
}

############################

rule_bsFs <- rule_closure(n, method = binSeg_fixedSteps)
rule_flFs <- rule_closure(n, method = fLasso_fixedSteps)
criterion <- function(x, vec){x}

bsFs_0JumpPValue <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)
flFs_0JumpPValue <- simulationGenerator(rule_flFs, paramMat, criterion,
  trials, cores)

save.image(file = paste0("res/pvalue_noJump_", Sys.Date(), ".RData"))
quit(save = "no")
