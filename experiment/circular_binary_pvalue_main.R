#load

#setting
len_vec <- c(10, 50, 100, 500, 1000, 2500, 5000)
sigma_vec <- c(0.5, 1, 2)
trials <- 1000
quantile_vec <- seq(0, 100, length.out = 201)

#initializing
parameter_mat <- expand.grid(len_vec, sigma_vec)
res_mat <- matrix(0, ncol = nrow(parameter_mat), nrow = length(quantile_vec))
colnames(res_mat) <- apply(parameter_mat, 1, function(x){paste0("Len", x[1], ":Sig", x[2])})

#main_func
fun <- function(trial, vec){
  set.seed(trial)
  dat <- stats::rnorm(vec[1], vec[2])
  circular_bin(dat)
}

#run
for(i in 1:nrow(parameter_mat)){
  tmp <- foreach(trial=1:trials) %dopar% fun(trial, parameter_mat[i,])
  res_mat[,i] <- stats::quantile(tmp, probs = quantile_vec)
  
  save.image("../experiment/circular_binSeg_image.RData")
}

save(res_mat, file = "../experiment/circular_binSeg_table.RData")