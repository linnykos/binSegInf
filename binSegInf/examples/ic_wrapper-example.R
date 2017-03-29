set.seed(1)
y = c(nnrom(10),rnorm(10)+5,rnorm(10))
sigma = 1 
## obj = wildBinSeg_fixedSteps(y,10)
obj <- list(cp = c(10,20,25,13,17,5))##Make some fake changepoints
poly <- ic_wrapper(cp = obj$cp, y = y, consec=2, maxsteps=length(cp), sigma,
                   type = "bic", returntype = c("polyhedra"))
stoptime <- ic_wrapper(cp = obj$cp, y = y, consec=2, maxsteps=length(cp), sigma,
                   type = "bic", returntype = c("stoptime"))
print(poly)
print(stoptime)
