source("simulation_header.R")

trials <- 1000
n <- 100
cores <- 20

paramMat <- as.matrix(0)

rule_bsFs_closure <- function(n, gridsize = 100){
  function(vec){
    dat <- CpVector(n, 0, NA)
    y <- dat$data
    
    obj <- binSeg_fixedSteps(y, 1)
  
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
  
    res <- tryCatch({
      c(confidence_interval(y, poly, contrast, gridsize = gridsize), 0)
    }, warning = function(w){
      c(-Inf, Inf, -1)
    })
    
    c(abs(contrast %*% y), res[1:2], jumps(obj), res[3])
  }
}

############################

rule_bsFs <- rule_bsFs_closure(n, 50)
criterion <- function(x, vec){x}

bsFs_0JumpCI <- simulationGenerator(rule_bsFs, paramMat, criterion,
  trials, cores)

save.image(file = paste0("res/CI_noJump_bsFs_", Sys.Date(), ".RData"))
quit(save = "no")
