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

    fun <- function(y){set.seed(y); criterion(rule(paramMat[x,]), paramMat[x,])}
    if(is.na(cores)){
      sapply(1:trials, fun)
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
