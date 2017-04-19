load(file="./data/ng754-s11-sheet3.Rdata")

## Load data
chromenum = 10
chromename = paste0("Chromosome ", chromenum)
dat = alldat[[chromename]]
y = dat[,"Log2Ratio"]
y = y[!is.na(y)]

## Run algo + inference for 3 steps on it
sigma = get_sigma(y)
numSteps = 3
method <- binSeg_fixedSteps

obj <- method(y, numSteps)
poly <- polyhedra(obj)
p.bsfs = rep(NA,length(obj$cp))
names(p.bsfs) = obj$cp
contrast = list()
for(ii in 1:length(obj$cp)){
    contrast[[ii]] <- contrast_vector(obj, ii)
    p.bsfs[ii] <- poly.pval(y=y,
                            G=poly$ gamma,u=poly$u, v=contrast[[ii]],sigma=sigma, bits=100)$pv
}
p.bsfs = cbind(rep(isim,length(obj$cp)), obj$cp, p.bsfs)

## Save it
resid.cleanmn <- y - cleanmn
filename = "../data/coriell.Rdata"
y.orig<-y
save(y.orig,
     segment.means,
     segments,
     cleanmn,
     resid.cleanmn,
     newmn,
     newmn.short,
     std,
     file = filename)



## Helpers
get_sigma <- function(y){
  cps = sort(changepoints(wbs::sbs(y))$cpt.th[[1]])
  segment.inds = sapply(1:(length(cps)+1),
                      function(ii) (c(0,cps,length(y))[ii]+1):(c(0,cps,length(y))[ii+1]))
  mn = rep(NA,length(y))
  for(ind in segment.inds) mn[ind] <- mean(y[ind])
  return(sd(y-mn))
}
