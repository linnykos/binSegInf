source("simulation_header.R")

trials <- 500
n <- 100
cores <- 20

paramMat <- as.matrix(0)

rule_bsFs_closure <- function(n){
  function(void){
    y <- CpVector(n, 0, NA)$data
    
    obj <- binSeg_fixedSteps(y, 1)
  
    poly <- form_polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
    
    if(any(poly$gamma %*% y < poly$u)) stop()
  
    res <- pvalue(y, poly, contrast)
    c(res, get_jumps(obj))
  }
}

############################

rule_bsFs <- rule_bsFs_closure(n)
criterion <- function(x, vec){x}

bsFs_0JumpPValue <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)

save.image(file = paste0("res/pvalue_noJump_bsFs_", Sys.Date(), ".RData"))
quit(save = "no")
