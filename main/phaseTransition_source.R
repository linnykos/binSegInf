#library(devtools)
#install_github("linnylin92/binSegInf", ref = "kevin", subdir = "binSegInf")

library(binSegInf)
library(wbs)
library(genlasso)
library(foreach)
library(doMC)

simulationGenerator <- function(rule, paramMat, criterion, trials, 
 cores = NA){

  if(!is.na(cores)) registerDoMC(cores = cores)

  res <- lapply(1:nrow(paramMat), function(x){
    cat(paste0("Row ", x, " started!\n"))

    fun <- function(y){set.seed(y); criterion(rule(paramMat[x,]))}
    if(is.na(cores)){
      vec <- sapply(1:trials, fun)
      print(vec)
      vec
    } else {
      .adjustFormat(foreach(trial = 1:trials) %dopar% fun(trial))
    }
  })
  
  names(res) <- sapply(1:nrow(paramMat), function(x){
    paste0(paramMat[x,], collapse = "-")})
  
  res
}

.adjustFormat <- function(lis){
  len <- sapply(lis, length)
  if(length(unique(len)) != 1) return(lis)

  ncol <- length(unique(len))
  if(ncol == 1) return(as.numeric(unlist(lis)))
 
  vec <- as.numeric(unlist(lis))
  matrix(vec, ncol = ncol, byrow = T)
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

.oneJumpRuleClosure <- function(n, func){
  function(vec){
    func(CpVector(n, vec[1:2], vec[3])$data)
  }
}

