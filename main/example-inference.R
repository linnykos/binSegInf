## Setup
outputdir = "output"
library(RColorBrewer)
library(igraph)
library(stringdist)
source('funs.R')
## Source in the p-value functions from the selectiveInference package 
source('selectinf/selectiveInference/R/funs.inf.R')

# Replacement of all usages of pval.fl1d to poly.pval; 
pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
}

## Setting 1
## n = 4
## set.seed(0)
## y = c(rep(0,n/2),rep(5,n/2)) + rnorm(n,0,1)
## thresh = .1

## Setting 2
## n = 60
## sd = .5
## mn = c(rep(-3, n/4), rep(2, n/4), rep(-1, n/4), rep(1, n/4))
## set.seed(1)
## y = mn + rnorm(n,0,sd)

## Setting 3
n = 10 
set.seed(0)
sigma=1
lev1=0
lev2=0
y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)
thresh = 1
nsim=10000
pvals = rep(NA,nsim)
y = c(rep(lev1,n/2),rep(lev2,n/2)) + rnorm(n,0,sigma)

slist = elist = blist = Blist = zlist = Zlist = matrix(NA, nrow = n, ncol = 2^8)
    
## Do binary segmentation
binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n)

## Get Gamma matrix
ii = 1
Gammat = matrix(NA,nrow = nrow(trim(Blist))*n, ncol = n)
u = rep(NA, nrow(trim(Blist))*n)
nrowGammat = 0
for(my.j in 1:nrow(trim(Blist))){
    print(my.j)
    bsublist = Blist[my.j,]
    for(b in bsublist[!is.na(bsublist)]){
        cat("b is", b, fill=T)
        if(my.j==1){
            k=1
            my.se = c(1,n)
        } else {
            ## k = floor((which(b == bsublist)+1)/2)
            k = which(b==bsublist)
            my.se = c(slist[my.j-1,k], elist[my.j-1,k]) 
        }

        terminal = is.na(blist[my.j,k])
        my.halfspaces = halfspaces(s = my.se[1], e = my.se[2], b = b,
                                   thresh = thresh, n = n, y = y, terminal)
        newrows = my.halfspaces[["V"]]
        newconst = my.halfspaces[["u"]]
        if(dim(newrows)[1]==0) next
        newrowinds = nrowGammat + c(1:nrow(newrows))
        if(any(newrowinds > nrow(Gammat) )){
            Gammat = rbind(Gammat, matrix(NA,nrow=nrow(Gammat),ncol=n))
            u = c(u, rep(NA,length(u)))
        }
        Gammat[newrowinds,] = newrows
        u[newrowinds] = newconst 
        nrowGammat = nrowGammat + nrow(newrows)
        ii = ii + 1
    }
}

Gammat = trim(Gammat,"row")
u = trim(u)
head(Gammat)

bzlist = blist*zlist

## Do inference.
sort(collapse.prev(bzlist))

## Make conditional contrasts
test.b = 5
if(test.b %in% abs(bzlist)){
    v = rep(0,n)
    z = zlist[which(test.b == blist, arr.ind = T)]
    ends = c(0,sort(blist),n)
    ind = which(test.b == ends)
    my.se = ends[c(ind-1, ind+1)]
    left.b = (my.se[1]+1):(test.b)
    right.b = (test.b+1):(my.se[2])
    v[left.b] = -1/length(left.b)
    v[right.b] = 1/length(right.b)
    v = v*z
    pval    = pval.fl1d(y = y,
                        G = Gammat,
                        dik = v,
                        sigma,
                        u = u)
    pvals[isim] = pval
}
