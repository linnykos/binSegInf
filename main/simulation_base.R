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

  ncol <- unique(len)
  if(length(ncol) == 1 && ncol == 1) return(as.numeric(unlist(lis)))
  do.call(cbind, lis)
}
