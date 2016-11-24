## ##' function to do one-step wild binary segmentation
## ##' @param y data
## ##' @param k steps you would like to take (defaults to 1)
## ##' @param ni number of intervals
## ##' @param intervals optionally, give intervals 
## wild_binseg = function(y, ni, intvls=NULL){

##     ## Draw random intervals
##     n = length(y)
##     if(is.null(intvls)) { intervals = get.intervals(n, ni)} else{ intervals = intvls}
    

##     ## Initialize Gamma matrix
##     jk  = list()
##     G = matrix(NA, ncol = length(y), nrow = 1000)
##     G.jk = cplist(100+1) ## G.jk will store the row numbers corresponding
##     Gn = 0 

##     ## running maximizer
##     maximizers = data.frame(matrix(nrow=100,ncol=6))
##     names(maximizers) = c("maxloc","s","e","cusum","which.row","signs")

##     ## Calculate cusum statistics + collect halfspaces
##     for(ii in 1:ni){
##         s = intervals[ii,1]
##         e = intervals[ii,2]
##         yy = y[s:e]
##         cusums = getcusums(s = s ,
##                            e = e,
##                            y = c((if(s==1) c() else rep(0,s-1)),yy))
        
##         ## Characterize signs
##         signed.cusummat = (cusums$contrasts) * (cusums$signs)
##         signed.cusummat = cbind(signed.cusummat, matrix(0,nrow=nrow(signed.cusummat), ncol=n-e))
##         G[(Gn+1):(Gn+nrow(signed.cusummat)),] = signed.cusummat 


##         ## Store running maximizer
##         maximizers[ii,] = c(cusums$bmax,s,e, cusums$cusum, Gn+ cusums$bmax.cusums, cusums$signs[cusums$bmax.cusums]) 

##         ## Update row number
##         Gn = Gn+nrow(signed.cusummat)
##     }

##     ## Extract maximizer
##     mx = maximizers[which.max(maximizers[,"cusum"]),]
##     maximizing.row = mx[,"which.row"]
    
##     ## Characterize the maximizing row of /all/ rows
##     maximizing.cusummat =  - G[1:Gn,] +do.call(rbind, lapply(1:Gn, function(ii)G[maximizing.row,]))
##     stopifnot((maximizing.cusummat)%*%y>=0)
##     G[(Gn+1):(2*Gn),] = maximizing.cusummat
##     Gn = 2*Gn

##     return(list(G=G[1:Gn,], u = rep(0,Gn),y=y, cp = mx[,"maxloc"], dir = mx[,"signs"]))
## }


## ## Uniform p-values
## nsim=1000
## ps=rep(NA,nsim)
## for(isim in 1:nsim){
## sigma=1
## y = c(rnorm(5,0,sigma),rnorm(5,5,sigma))
## y = rnorm(10,0,sigma)
## n = length(y)
## output = wild_binseg(y,10)
## v = make.v.simple(b=output$cp,
##                   s=1,
##                   e=n,
##                   n,
##                   dir=output$dir)
## ps[isim] = poly.pval(y,output$G,output$u, v, sigma)$pv
## }
## hist(ps)


## ## For varying signal size
## nsim=100
## deltas = c(0,2,4)
## ps.rwbs = ps.wbs=ps.bs = matrix(NA,nrow=nsim, ncol = length(deltas))
## sigma=1
## n=10
## for(idelta in 1:length(deltas)){
##     print(idelta)
##     delta = deltas[idelta]
##     for(isim in 1:nsim){
##         print(isim)
        
##         ## Generate data
##         y = c(rnorm(n/2,0,sigma),rnorm(n/2,delta,sigma))



##         ## Wild binseg
##         output = wild_binseg(y,10)
##         if(output$cp==n/2){
##             v = make.v.simple(b=output$cp, s=1, e=n, n, dir=output$dir)
##             ps.wbs[isim,idelta] = poly.pval(y,output$G,output$u, v, sigma)$pv
##         }

##         ## Binseg
##         a = binseg.by.size(y, numsteps=1, verbose=FALSE)
##         if(a$B[1]==n/2){
##             v = make.v(a$B[1],a$B,a$Z,n)
##             ps.bs[isim,idelta]=poly.pval(y,a$G,a$u,v,sigma)$pv
##         }


##         ## Randomized wild binseg
##         do.randomized.inference = function(y,iwsim){
##             ps = rep(NA,iwsim)
##             for(iw in 1:iwsim){
##                 w = get.intervals(n,ni)
##                 ## y = c(rnorm(n/2,0,sigma),rnorm(n/2,delta,sigma))
##                 output = wild_binseg(y,10,w)
##                 if(output$cp==n/2){
##                     v = make.v.simple(b=output$cp, s=1, e=n, n, dir=output$dir)
##                     ps[iw] = poly.pval(y,output$G,output$u, v, sigma)$pv
##                 }
##             }
##             return(mean(ps[!is.na(ps)]))
##         }
##         p.rwbs = do.randomized.inference(y,100)
##         ps.rwbs[isim,idelta] = p.rwbs
##     }
## }



## ## Make qqplots
## pcols = brewer.pal(100,"Set3")
## ii=3
## make.qqplot.background()
## make.qq.line(p=ps.bs[,ii],pch=16,pcol=pcols[3])
## make.qq.line(p=ps.wbs[,ii],pch=16,pcol=pcols[4])
## make.qq.line(p=ps.rwbs[,ii],pch=16,pcol=pcols[5])
## title(bquote(delta==.(deltas[ii])))
## legend("bottomright",legend=c("BS","WBS","RWBS"),pch=c(16,16,16),col=pcols[3:5])

