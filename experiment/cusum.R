trials <- 10000
n <- 100

mat <- matrix(0, trials, n-1)
for(trial in 1:trials){
  set.seed(10*trial)
  y <- rnorm(n)
  
  breakpoint <- seq(from = 1, to = n - 1, by = 1)
  #mat[trial,] <- sapply(breakpoint, binSegInf:::.cusum, y = y, start = 1, end = n)
  mat[trial,] <- sapply(breakpoint, function(x){
    sqrt(1/(1/x + 1/(n-x)))*(mean(y[(x+1):n]) - mean(y[1:x]))
  })
  
  if(trial %% floor(trials/10) == 0) cat('*')
}

mean.vec <- apply(mat, 2, mean)
var.vec <- apply(mat, 2, var)

par(mfrow = c(1,2))
plot(mean.vec)
plot(var.vec)

break.vec <- apply(mat, 1, function(x){which.max(abs(x))})
hist(break.vec, breaks = 20)

par(mfrow = c(1,3))
hist(mat[,1], breaks = 20, col = "gray")
hist(mat[,50], breaks = 20, col = rgb(1,0,0,0.5), add = T)

plot(mat[,1], mat[,50], pch = 16)
lines(x=c(-10,10), y = c(-10,10), col = "red")

idx <- which(break.vec == 1)
plot(mat[idx,1], mat[idx,50], pch = 16)
lines(x=c(-10,10), y = c(-10,10), col = "red")


#############################################

#BINARY SEG
cov.mat <- matrix(0, n-1, n-1)
for(i in 1:(n-1)){
  for(j in i:(n-1)){
    w1 <- sqrt(1/(1/i + 1/(n-i)))
    w2 <- sqrt(1/(1/j + 1/(n-j)))
    cov.mat[i,j] <- w1*w2*(1/j - (j-i)/((n-i)*j) + 1/(n-i))
    cov.mat[j,i] <- cov.mat[i,j]
  }
}

clockwise90 = function(a) { t(a[nrow(a):1,]) }
image(clockwise90(cov.mat))

break.vec.sim <- numeric(trials)
for(trial in 1:trials){
  set.seed(10*trial)
  vec <- MASS::mvrnorm(1, rep(0, n-1), cov.mat)
  break.vec.sim[trial] <- which.max(abs(vec))
  
  if(trial %% floor(trials/10) == 0) cat('*')
}
hist(break.vec.sim)

####################################

#FUSED LASSO
cov.mat <- matrix(0, n-1, n-1)
for(i in 1:(n-1)){
  for(j in i:(n-1)){
    cov.mat[i,j] <- i*(n-j)/n
    cov.mat[j,i] <- cov.mat[i,j]
  }
}

image(clockwise90(cov.mat))

break.vec.sim <- numeric(trials)
for(trial in 1:trials){
  set.seed(10*trial)
  vec <- MASS::mvrnorm(1, rep(0, n-1), cov.mat)
  break.vec.sim[trial] <- which.max(abs(vec))
  
  if(trial %% floor(trials/10) == 0) cat('*')
}
hist(break.vec.sim)

#####################
y <- rnorm(n)
D <- binSegInf:::.form_Dmatrix(n)

trials <- 1000
mat <- matrix(0, trials, n-1)
for(trial in 1:trials){
  set.seed(10*trial)
  y <- rnorm(n)
  
  mat[trial,] <- binSegInf:::.compute_fused_numerator(D, 1:(n-1), y)
  
  if(trial %% floor(trials/10) == 0) cat('*')
}

mean.vec <- apply(mat, 2, mean)
var.vec <- apply(mat, 2, var)

par(mfrow = c(1,2))
plot(mean.vec)
plot(var.vec)

