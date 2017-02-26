binSeg_cov <- function(n){
  mat <- matrix(0, n-1, n-1)
  for(i in 1:(n-1)){
    for(j in i:(n-1)){
      weight <- sqrt(1/((1/i + 1/(n-i))*(1/j+1/(n-j))))
      mat[i,j] <- weight *(1/j - (j-i)/(j*(n-i)) + 1/(n-i))
      mat[j,i] <- mat[i,j]
    }
  }
  
  mat
}

fLasso_cov <- function(n){
  mat <- matrix(0, n-1, n-1)
  for(i in 1:(n-1)){
    for(j in i:(n-1)){
      mat[i,j] <- i*(n-j)/n
      mat[j,i] <- mat[i,j]
    }
  }
  
  mat
}

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

###############################
#plot
n <- 50
bsCov <- binSeg_cov(n); flCov <- fLasso_cov(n)

par(mfrow = c(1,2))
image(.rotate(bsCov), asp = T); image(.rotate(flCov), asp = T)

##############################

bsEig <- eigen(bsCov); flEig <- eigen(flCov)
plot(bsEig$values, pch = 16); plot(flEig$values, pch = 16)

bsVec <- bsEig$vectors[,1]; flVec <- flEig$vectors[,1]
plot(bsVec, pch = 16); plot(flVec, pch = 16)

.l2norm <- function(vec){
  sqrt(sum((vec)^2))
}

bsVec %*% flVec / (.l2norm(bsVec) * .l2norm(flVec))
acos(bsVec %*% flVec / (.l2norm(bsVec) * .l2norm(flVec)))

################################

#try a simple case of n = 4
n <- 4; trials <- 1000
bsCov <- binSeg_cov(n); flCov <- fLasso_cov(n)
set.seed(10); bsRes <- apply(abs(MASS::mvrnorm(trials, rep(0, n-1), bsCov)), 1, which.max)
set.seed(10); flRes <- apply(abs(MASS::mvrnorm(trials, rep(0, n-1), flCov)), 1, which.max)

table(bsRes); table(flRes)
plot(table(bsRes)); plot(table(flRes))

###############################

library(rgl)
rgl::plot3d(rgl::ellipse3d(bsCov))
rgl::plot3d(rgl::ellipse3d(flCov))


########################

#would the results change if we zeroed out all the other eigenvalues?
n <- 50; trials <- 1000; trunc <- 3
bsCov <- binSeg_cov(n); flCov <- fLasso_cov(n)
bsEig <- eigen(bsCov); flEig <- eigen(flCov)
bsEig$values[-c(1:trunc)] <- 0; flEig$values[-c(1:trunc)] <- 0
bsCov <- bsEig$vectors%*%diag(bsEig$values)%*%t(bsEig$vectors)
flCov <- flEig$vectors%*%diag(flEig$values)%*%t(flEig$vectors)

set.seed(10); bsRes <- apply(abs(MASS::mvrnorm(trials, rep(0, n-1), bsCov)), 1, which.max)
set.seed(10); flRes <- apply(abs(MASS::mvrnorm(trials, rep(0, n-1), flCov)), 1, which.max)

par(mfrow = c(1,2))
plot(table(bsRes)); plot(table(flRes))

