library(binSegInf)
library(wbs)
library(genlasso)

simulationGenerator <- function(rule, paramMat, criterion, trials){
  res <- lapply(1:nrow(paramMat), function(x){
    sapply(1:trials, function(y){set.seed(y); criterion(rule(x))})
  })
  
  names(res) <- sapply(1:nrow(paramMat), function(x){
    paste0(paramMat[x,], collapse = "-")})
  
  res
}

.noJumpRuleClosure <- function(n, func){
  function(void){func(CpVector(n, 0, NA)$data)}
}

.extractJumps <- function(obj, ...) UseMethod(".extractJumps")

.extractJumpClosure <- function(numJumps){
  function(res){.extractJumps(res, numJumps)}
}

.extractJumps.fusedlasso <- function(obj, numJumps, ...){
  enumerateJumps(coef(obj, df = numJumps + 1)$beta)
}

.extractJumps.sbs <- function(obj, numJumps, ...){
  order(obj$res[,"min.th"], decreasing = TRUE)[1:numJumps]
}

#####################################

noJumpRuleFl <- .noJumpRuleClosure(100, fusedlasso1d)
noJumpRuleBs <- .noJumpRuleClosure(100, sbs)
noJumpCriterion <- .extractJumpClosure(1)
paramMat <- matrix(0, 1, 1)
trials <- 10

resFl <- simulationGenerator(noJumpRuleFl, paramMat, noJumpCriterion, trials)
resBs <- simulationGenerator(noJumpRuleBs, paramMat, noJumpCriterion, trials)



