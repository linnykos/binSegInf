## Synopsis: The idea is to add bootstrapped residuals, then see the 1.
## conditional/unconditional powers of the tests done at the viscinities of the
## true guys, and 2. uniformity (or lack thereof) in the tests conducted
## regarding null locations.

# Data directory
datadir = "../data"
filename = "coriell05296.Rdata"
load(file=file.path(datadir,filename))
source(file=file.path("../main/artificial/artif-helpers.R"))


## Factor
## fac = 2 ## 1, 1.5, 2, 2.5, 3
whichfac = 5
fac = c(1, 1.5, 2, 2.5, 3)[whichfac]
nsim = c(200,300,400,500,600)[whichfac]
## nsim = 200
bits = 3000
mc.cores = 7
results = list()
## for(isim in 1:nsim){
start.time = Sys.time()
results = mclapply(1:nsim, function(isim){

    printprogress(isim, nsim,
                  lapsetime = round(difftime(Sys.time(), start.time,
                                             units = "secs"), 2))
    ## Run a single result
    myresult = onesim_rbs(y.orig, bits=bits, fac=fac, verbose=FALSE)
}, mc.cores=mc.cores)




## Frequency of captures
freqs.wid5=freqs.wid1=c()
for(ifac in c(1,3,4,5)){
    fac = c(1, 1.5, 2, 2.5, 3)[ifac]
    facstring = paste0(unlist(strsplit(toString(fac), split='.', fixed=TRUE)), collapse="")
    filename = paste0("artif-rbs-fac-", facstring, ".Rdata")
    load(file=file.path(outputdir, filename))

    ## extract p-values
    results.proper = results[sapply(results, length)!=1]
    pvlist = sapply(results.proper, function(a)a$pvs)

    ## Get all locations
    all.locs = abs(as.numeric(unlist(sapply(pvlist, names))))
    newmn2 = newmn[-(1:200)]
    true.locs = which(newmn2[2:length(newmn2)] - newmn2[1:(length(newmn2)-1)]!=0)


    ## Get frequencies
    wid=1
    true.locs.vic = sapply(true.locs, function(myloc)myloc + seq(from=-wid,to=wid))
    freqs.wid1[ifac] = (sum(all.locs %in% true.locs.vic)/length(all.locs))
    wid=5
    true.locs.vic = sapply(true.locs, function(myloc)myloc + seq(from=-wid,to=wid))
    freqs.wid5[ifac] = (sum(all.locs %in% true.locs.vic)/length(all.locs))

    ## Separate true locatiosn and false ones.
    nonnull.inds = which(abs(as.numeric(names(unlist(pvlist)))) %in% as.numeric(true.locs.vic))
    null.pvs = unlist(pvlist)[-nonnull.inds]
    qqunif(null.pvs)
}

## Plot frequency of captures
facs = c(1, 1.5, 2, 2.5, 3)
pdf(file=file.path(outputdir, "freq.captured"))
plot(freqs.wid5~facs, cex=2, pch=16, type='p' , ylim=c(0,1))
title(main="How frequently are locations captured \n in +-5 (black), +-1 (red) neighborhood, overall")
lines(freqs.wid1~facs, cex=2, pch=16, type='p', ylim=c(0,1), col='red')
graphics.off()

## Data plot
plot(y)
lines(density(all.true.locs)$y*100~density(all.true.locs)$x, lwd=3)
lines(newmn[-(1:200)], lwd=3, col='red')
abline(v=true.locs)


## Try forming /actual contrasts/ to see this.
null.pvs.list = list()
numtests.list = list()
null.rejects.list = list()
null.fwer.list = list()
for(ifac in c(1,3,4,5)){
    fac = c(1, 1.5, 2, 2.5, 3)[ifac]
    facstring = paste0(unlist(strsplit(toString(fac), split='.', fixed=TRUE)), collapse="")
    filename = paste0("artif-rbs-fac-", facstring, ".Rdata")
    load(file=file.path(outputdir, filename))
    results.proper = results[sapply(results, length)!=1]

    null.pvs.list = sapply(results.proper, function(myproperresult){
        cp = (myproperresult$locs.retained)
        cp.sign = sign(cp)
        cp = abs(cp)
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(newmn[-(1:200)]))
        means = sapply(vlist, function(v)v %*% newmn[-(1:200)])
        null.cp = cp[which(means==0)]
        pvs = myproperresult$pvs
        pvs = pvs[which(abs(as.numeric(names(pvs))) %in% null.cp)]
        return(pvs)
    })


    ## The total number of tests
    numtests = sapply(results.proper, function(myproperresult){
        cp = (myproperresult$locs.retained)
        numtests = length(cp)
    })


    ## Calculate things
    null.rejects = sapply(null.pvs.list, function(null.pvs)null.pvs<0.05)
    null.fwer = Map(function(null.pvs, numtest){any(null.pvs<0.05/numtest)},
                    null.pvs.list, numtests)

    ## Store them
    null.pvs.list[[ifac]] = unlist(null.pvs)
    null.rejects.list[[ifac]] = unlist(null.rejects)
    null.fwer.list[[ifac]] = unlist(null.fwer)
    numtests.list[[ifac]] = numtests
}

## Calculate errors
type.I.errors = sapply(c(1,3,4,5), function(ifac){
    sum(null.rejects.list[[ifac]])/length(null.rejects.list[[ifac]]) })
fwer = sapply(c(1,3,4,5), function(ifac){
    sum(null.fwer.list[[ifac]])/length(null.fwer.list[[ifac]])
})

## Make plot
jpeg(file=file.path(outputdir, "artif-rbs-errors.jpg"), width=500, height=500)
plot(type.I.errors~c(1,3,4,5), ylim = c(0,1), pch=16, cex=1)
points(fwer~c(1,3,4,5), ylim = c(0,1), pch=16, cex=1, col='red')
title(main="Type I errors (black) \n and FWERs (red)")
graphics.off()
## Not sure if type I errors

sum(null.rejects.list[[1]])/sum(total.numtests.list[[1]])

## Individual null p-values
## pdf(file=file.path(outputdir,"artif-rbs-nullpvals-by-fac.pdf"), width=7,height=7)
jpeg(file=file.path(outputdir,"artif-rbs-nullpvals-by-fac.jpg"), width=700,height=700)
par(mfrow=c(2,2))
for(ifac in c(1,3,4,5)){
    qqunif(null.pvs.list[[ifac]])
    title(paste0("Null p-values for fac=",facs[ifac]))
}
graphics.off()

## Overall
## pdf(file=file.path(outputdir,"artif-rbs-nullpvals-combined.pdf"), width=7,height=7)
jpeg(file=file.path(outputdir,"artif-rbs-nullpvals-combined.jpg"), width=500,height=500)
qqunif(unlist(null.pvs.list))
title(main="Overall null p-values")
graphics.off()

sum(unlist(null.pvs.list)<0.05)/length(unlist(null.pvs.list))


## Noise distribution
jpeg(file=file.path(outputdir,"artif-rbs-noise1.jpg"), width=500,height=500)
par(mfrow=c(2,1))
plot(y.orig[1:200], main = "first 200 points")
qqnorm(y.orig[1:200])
qqline(y.orig[1:200])
graphics.off()


## SDs over time
sds = sapply(seq(from=1,to=900, by=10),function(istart){
    sd(y.orig[istart:(istart+200)])
})
jpeg(file=file.path(outputdir,"artif-rbs-noise2.jpg"), width=500,height=500)
par(mfrow=c(2,1))
plot(y.orig)
abline(v=c(0,200), lwd=2, col='blue')
abline(v=c(900,1100), lwd=2, col='blue')
plot(sds)
graphics.off()



## bootstrap distribution in first 1100 points
y[1:1100]


