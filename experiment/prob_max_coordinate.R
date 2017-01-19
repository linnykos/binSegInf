probability_max_coordinate <- function(mean, Sigma, i){
  stopifnot(i < length(mean), is.matrix(Sigma), all(length(mean) == dim(Sigma)))
  
  #compute the n-2 different means
  mean_dif <- mean[i] - mean[-i]
  
  #compute the n-2 different covariances
  cov_dif <- Sigma[i,i] + diag(Sigma)[-i] - 2*Sigma[i,-i]
  
  #add their log of the normal cdfs together
  log_prob <- sapply(1:length(mean_dif), function(x){
    log(1-pnorm(0, mean = mean_dif[x], cov_dif[x]))
  })
  
  #exponeniate
  exp(sum(log_prob))
}

probability_max_coordinate(rep(0, n-1), cov.mat1, 50)
