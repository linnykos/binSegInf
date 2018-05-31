## sim = function(n=100, nrep=5000, err.fun=rnorm, mu=rep(0,n), seed=NULL) {
fl_onestep <- function(y){
    n=length(y)
    
    G = matrix(0,2*n-4,n)
    M = A = matrix(0,n-1,n)
    fac = rep(NA, n-1)
    for (i in 1:(n-1)) {
      M[i,] = c(rep(-1/i,i),rep(1/(n-i),n-i))
      A[i,] = i*(n-i)/n * M[i,]
      fac[i] = i*(n-i)/n
    }
  
    ## Fused lasso
    j1 = which.max(abs(A %*% y))
    s1 = sign(sum(A[j1,]*y))
    G[1:(n-2),] = t(s1*A[j1,] - t(A[-j1,]))
    G[(n-1):(2*n-4),] = t(s1*A[j1,] + t(A[-j1,]))
    obj = pval(y, G, s1*M[j1,])
    j.fl = j1
    p.fl = obj$p
    vlo.fl = obj$vlo
    vup.fl = obj$vup

    return(list(cp=j1, cp.sign=s1, pv=p.fl, A=A, M=M, fac=fac))
}

## Does this match the (detection and inference) first-step fused lasso? Yes!
set.seed(0)
D = makeDmat(n,order=0)
nrep = 100
for(isim in 1:nrep){
    printprogress(isim,nrep)
    y = rnorm(n,0,1)
    out.new = fl_onestep(y)
    p.new = out.new$pv

    ## Run the old way
    out = dualpathSvd2(y, D=D, maxsteps=1)
    G = polyhedra(out)
    i = out$cp
    v = c(rep(-1/i,i),rep(1/(n-i),n-i))
    v = v * out$cp.sign
    p = pval(y, G$gamma, v)$p
    
    ## Match several things
    stopifnot(out.new$cp * out.new$cp.sign == out$cp * out$cp.sign) 
    stopifnot(all.equal(p, p.new))
}





## by 12:15p, have a script of the coefficients ready.



## Step 2
Aold = out$coeflist[[1]]

plot(A%*%y)

plot(Aold[10,])
lines(Aold.binseg[10,],col='red')
abline(h=0)

plot((out$coeflist[[1]]/out$base.coeflist[[1]])[,1])
plot(abs(out$coeflist[[2]]/out$base.coeflist[[2]])[,1])
plot(abs(out$coeflist[[2]]/out$base.coeflist[[2]])[,60])


G = polyhedra(out)
i = out$cp
v = c(rep(-1/i,i),rep(1/(n-i),n-i))
v = v * out$cp.sign
p = pval(y, G$gamma, v)$p


## 
