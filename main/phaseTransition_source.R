.noJumpRuleClosure <- function(n, func){
  function(void){func(CpVector(n, 0, NA)$data)}
}

.extractJumps <- function(obj, ...) UseMethod(".extractJumps")

.extractJumpClosure <- function(numJumps){
  function(res){.extractJumps(res, numJumps)}
}

.extractJumps.fusedlasso <- function(obj, numJumps, ...){
  enumerateJumps(coef(obj, df = numJumps + 1)$beta)[1:numJumps]
}

.extractJumps.sbs <- function(obj, numJumps, ...){
  order(obj$res[,"min.th"], decreasing = TRUE)[1:numJumps]
}

.oneJumpRuleClosure <- function(n, func){
  function(vec){
    func(CpVector(n, vec[1:2], vec[3])$data)
  }
}

