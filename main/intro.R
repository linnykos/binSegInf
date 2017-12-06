## Synopsis: introductory example for paper
load(file="../data/ng754-s11-sheet3.Rdata")
library(xtable)

## Load data
chromenum = 10
chromename = paste0("Chromosome ", chromenum)
dat = alldat[[chromename]]
y = dat[,"Log2Ratio"]
y = y[!is.na(y)]

## Run algo + inference for 3 steps on it
sigma = get_sigma(y) ## Harvested from a reasonable binary segmentation.
numSteps = 4
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

## Naive inference
p.naive = sapply(1:length(obj$cp), function(ii) ztest(contrast[[ii]], y,sigma,0.05))
names(p.naive) = names(p.bsfs)

xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-0.75,1.5)
mar = c(4.5,4.5,0.5,0.5)
cp = sort(obj$cp)
beta0 <- get_piecewise_mean(y, cp)
beta1 <- get_piecewise_mean(y, c(53,95))
Letters = toupper(letters[1:length(cp)])

# Plot, data and mean
## pdf(file=file.path("$DB/research/binSegInf/writing/figures", "intro1.pdf"),width=w,height=h)
pdf(file=file.path("../main/figures", "intro.pdf"),width=w,height=h)
par(mar=mar)
plot(y, ylim = ylim, axes=F, xlab = xlab, ylab = "", pch=pch, col=pcol);
lines(beta0, col='blue', lwd=lwd)
lines(beta1, col='red', lwd=lwd)
text(Letters, x = cp+2, y = rep(1.2,1.2))
abline(v = cp, col = "black", lty=2)
axis(1); axis(2)
legend("bottomright", pch=c(pch,NA,NA,NA), lty=c(NA,1,1,2),lwd=c(NA,lwd,lwd,1),
       col = c(pcol,"blue","red","black"),legend=c("aCGH Data","Estimate","Estimate for simulation","Changepoint"), bg="white")
graphics.off()

# Generate table of p-values
mytable = cbind(p.naive,p.bsfs)
mytable = signif(mytable,3)
mytable = mytable[order(obj$cp),]
mytable = cbind(rownames(mytable),mytable)
rownames(mytable) = Letters
colnames(mytable) = c("Location","Naive p-values","TG p-values")
xtable(mytable, digits=c(1,0,3,3), align="|r|r|r|r|")



## Bootstrap residuals and visualize
cleanmn <- piecewise_mean(y, c(53,94))
std = sd(y-cleanmn)
par(mfrow=c(3,3))
plot(y)
lines(cleanmn,col='red')
for(ii in 1:8){
plot(cleanmn + bootstrap_sample(y-cleanmn))
lines(cleanmn,col='red')
}

##  Bootstrap residuals
onesim_intro <- function(isim, bootstrap=TRUE){
  print(isim)
  numSteps = 3
  method <- binSeg_fixedSteps
  seed = NULL##isim
  std = sd(y-cleanmn)
  if(bootstrap) y = (cleanmn + bootstrap_sample(y-cleanmn,seed=seed))
  if(!bootstrap) y = (cleanmn + rnorm(length(cleanmn),0,std))
  obj <- method(y, numSteps)
  poly <- polyhedra(obj)
  p.bsfs = rep(NA,length(obj$cp))
  names(p.bsfs) = obj$cp
  contrast = list()
  for(ii in 1:length(obj$cp)){
    contrast[[ii]] <- contrast_vector(obj, ii)
    p.bsfs[ii] <- poly.pval(y=y,
                            G=poly$ gamma,u=poly$u, v=contrast[[ii]],sigma=std, bits=100)$pv
  }
  p.naive = sapply(1:length(obj$cp), function(ii){ztest(contrast[[ii]], y,std)})
  names(p.naive) = names(p.bsfs)
  return(list(p.bsfs = p.bsfs, p.naive = p.naive))
}


p.bsfs = c()
nsim=1000
manysimresult = lapply(1:nsim, onesim_intro)
manysimresult.Gaussian= lapply(1:nsim, onesim_intro, FALSE)


## Extract pvalue list
plist.bsfs = lapply(manysimresult, function(a)(a$p.bsfs))
list.naive = lapply(manysimresult, function(a)(a$p.naive))

## Extract spurious p-vals
get_spurious <- function(my.pmat,cp,slack){
n = ncol(my.pmat)
  a = lapply(cp, function(pt){
    visc <-  pt + (-slack):slack
    tested.at.all = apply(my.pmat, 1, function(myrow){
      tested = any(!is.na(myrow[visc]))
    })
    return(tested.at.all)
  })
interesting.rows = apply(do.call(cbind,a),1, all)
all.visc = unique(do.call(c,lapply(cp, function(pt){pt+(-slack):(slack)})))
all.outside.visc =(1:n)[-all.visc]
spurious.pmat = (my.pmat[interesting.rows, all.outside.visc])
spurious.pvals = as.numeric(spurious.pmat)
spurious.pvals = spurious.pvals[!is.na(spurious.pvals)]
}

## Viscinity
slack=5
my.pmat = reformat(plist.bsfs,n)
pv1 = get_condit_pvals(my.pmat, loc = 53 + (-slack):slack)
pv2 = get_condit_pvals(my.pmat, loc = 95 + (-slack):slack)
pv3 = get_spurious(my.pmat, c(53,95), slack)


## Reformat naive p-values that missed
my.pmat = reformat(plist.naive,n)
pv3v2 = get_spurious(my.pmat, c(53,95), slack)## get_condit_pvals(my.pmat, (1:ncol(my.pmat))[-visc] )

## Make plot
w=5;h=5
filename = "bootstrap-intro.pdf"
pdf(file.path("../main/figures",filename), width=w,height=h)
xlab="Observed"
ylab="Observed"
xlim=c(0,1)
ylim=c(0,1)
mar = c(4.5,4.5,0.5,0.5)
cols = RColorBrewer::brewer.pal(4,"Set3")
par(mar=mar)
plot(NA,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab)
abline(0,1)
axis(2);axis(1)
Map(function(pv,pch,col){qqunif_add(pp=pv, pch=pch, col=col)}, list(pv1,pv2,pv3,pv3v2), rep(pch,4), cols)
legend("bottomright",
       col=cols,
       pch = rep(16,4),
       lwd = rep(1,4),
       ## legend=c("A","B","C","D", "other", "naive other"))
       legend=c("A","C","other", "naive other"))
graphics.off()


## What if I compare this to having drawn Normal noise, with the same noise level? Same deal!
plist.bsfs.Gaussian  = lapply(manysimresult.Gaussian, function(a)(a$p.bsfs))
plist.naive = lapply(manysimresult.Gaussian, function(a)(a$p.naive))
