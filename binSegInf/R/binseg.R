#' The main binary segmentation function. When \code{thresh} is set to zero,
#' then this can be used to collect the unbalanced haar wavelet basis.
#' @param s Starting index, in the vector-valued data. Must be an integer larger
#'     than or equal to 1, and strictly smaller than \code{e}.
#' @param e Ending index, in the vector-valued data. Must be an integer smaller
#'     than or equal to \code{n}, and strictly larger than \code{s}.
#' @param j The depth of the recursion on hand.
#' @param k The indexing of the node location, from left to right, in the
#'     \emph{complete} binary tree.
#' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
#'     for the recursion.
#' @param y The original data.
#' @param n The length of the data \code{y}.
#' 
#' @examples
#' n = 60
#' sd = .5
#' mn = c(rep(-3, n/4), rep(2, n/4), rep(-1, n/4), rep(1, n/4))
#' y = mn + rnorm(n,0,sd)
#' ## blist are the breaks, and zlist are the associated signs.
#' blist = zlist = matrix(NA, nrow = n, ncol = 2^8)
#' thresh = .5
#' ##binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n) ## check() doesn't like this; not sure why.
binseg = function(s, e, j, k, thresh, y, n, verbose=F){
    if(verbose){
        cat("s,e,j,k",fill=T)
        cat(s,e,j,k,fill=T)
    }
    if(e-s<1){
       if(verbose) cat("terminated because e-s<1", fill=T)
       slist[j,k] <<- s
       elist[j,k] <<- e
       return() 
    } else {
        all.bs = (s:(e-1))
        df = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
        df = c(rep(NA,s-1),df, rep(NA,n-s))
        sn = sign(df) 

        ## Obtain breakpoint and its sign
        b = which.max(abs(df))
        z = sn[b]
        if(verbose) cat("b is", b, fill=T)

        ## Check threshold exceedance, then store
        if(abs(df[b]) < thresh){
            Blist[j+1,k] <<- b
            Zlist[j+1,k] <<- z
            slist[j,k] <<- s
            elist[j,k] <<- e
            if(verbose) cat("terminated because biggest gap was",abs(df[b]),fill=T)
            return()
        } else { 
            if(verbose) cat(df[b], fill=T)
            blist[j+1,k] <<- b
            Blist[j+1,k] <<- b
            zlist[j+1,k] <<- z
            Zlist[j+1,k] <<- z
            slist[j,k] <<- s
            elist[j,k] <<- e
        }
                    
        ## Recurse
        binseg(s, b, j+1, 2*k-1, thresh, y, n)
        binseg(b+1, e, j+1, 2*k, thresh, y, n)
    }
}





##' Function to carry out fixed-number-of-steps binary segmentation.
##' @param y numeric vector, data
##' @param numsteps desired number of changepoints
##' @param verbose set to true for algorithm run details. 
#### ' @import Matrix

binseg.by.size = function(y,numsteps,verbose=FALSE){

## Basic checks
    if(numsteps > length(y)-1) stop(paste("You should ask for less than", length(y), "steps!"))
    
    ## Initialize things
    require(Matrix)   
    B = Z = rep(NA,length(y)) 
    Bcurr = Zcurr = Scurr = Ecurr =
        Matrix(0, ncol=2^(numsteps+1), nrow = numsteps+1, sparse=TRUE)
    S = E = A = Tt =  
        Tcurr = Acurr =
            lapply(1:length(y),function(i) c())
    Scurr[1,1] = 1
    Ecurr[1,1] = length(y)
    Tcurr[[1]] = c(1,1)
    Acurr[[1]] = c()
    jk  = list()
    zetas = rep(NA,length(y))
    G = matrix(NA, ncol = length(y), nrow = 2*length(y)*numsteps)
    Gn = 0 
    
    ## Main loop
    for(mystep in 1:numsteps){
        if(verbose) cat("At step", mystep, " ") ## Get all candidate changepoints
        curr.max=-Inf
        Gn.beginning.of.step = Gn
        for(ii in 1:sum(!sapply(Tcurr, is.null))){ 
            Gn.beginning.of.this.node = Gn
            j = Tcurr[[ii]][1]
            k = Tcurr[[ii]][2]
            if(Ecurr[j,k]-Scurr[j,k]<=1) next
            cusums = getcusums(s = Scurr[j,k],
                               e = Ecurr[j,k],
                               y = y)
            
            ## Characterize signs
            signed.cusummat = (cusums$contrasts) * (cusums$signs)
            G[(Gn+1):(Gn+nrow(signed.cusummat)),] = signed.cusummat 
            Gn = Gn+nrow(signed.cusummat)
            
            ## Find cusum maximizer
            Bcurr[j, k] = cusums$bmax
            Zcurr[j, k] = cusums$signs[cusums$bmax.cusums]
            breaking.cusum = cusums$cusum
            
            ## Keep running maximum
            if(curr.max <= breaking.cusum){
                curr.which.max = c(j,k) 
                curr.max = breaking.cusum
                curr.max.signed.row = signed.cusummat[cusums$bmax.cusums,]
                curr.max.signed.rownum = Gn.beginning.of.this.node + cusums$bmax.cusums
            }
        }
        
        ## Record maximizer row and rownum
        max.signed.row = curr.max.signed.row
        max.signed.rownum = curr.max.signed.rownum
        
        ## Characterize cusum-maximizer
        this.step.rows=(Gn.beginning.of.step+1):Gn
        this.step.rows = this.step.rows[this.step.rows!=max.signed.rownum]
        comparison.cusummat = t(apply(G[this.step.rows,,drop=FALSE], 1, function(myrow){
            curr.max.signed.row - myrow }))
        G[(Gn+1):(Gn+nrow(comparison.cusummat)),] = comparison.cusummat
        Gn = Gn+nrow(comparison.cusummat)

        ## Check characterization (to be continued)


        
        ## Record knot as CUSUM maximizer
        zetas[mystep] = curr.max
        jmax = curr.which.max[1]
        kmax = curr.which.max[2]
        
        ## Update terminal and active node set
        which.duplicate = which(sapply(Tcurr, function(myjk){all.equal(myjk, c(jmax, kmax))==TRUE}))
        Tcurr[[which.duplicate]] <- c(jmax + 1, 2*kmax - 1)
        Tcurr[[mystep+1]] <- c(jmax + 1, 2*kmax)
        Acurr[[mystep+1]] <- c(jmax,kmax)
        
        ## Update Scurr and Ecurr for the /new/ nodes
        Scurr[jmax+1,2*kmax-1] = Scurr[jmax,kmax]
        Ecurr[jmax+1,2*kmax-1] = Bcurr[jmax,kmax]
        
        Scurr[jmax+1,2*kmax] = Bcurr[jmax,kmax]+1
        Ecurr[jmax+1,2*kmax] = Ecurr[jmax,kmax]
        
        ## Take snapshot
        S[[mystep]] = trim(Scurr)
        E[[mystep]] = trim(Ecurr)
        A[[mystep]] = trim(Acurr)
        Tt[[mystep]] = trim(Tcurr)
        B[mystep] = Bcurr[jmax,kmax]
        Z[mystep] = Zcurr[jmax,kmax]
        jk[[mystep]] = c(jmax,kmax) 
        
        if(verbose) cat("breakpoint", Bcurr[jmax,kmax], "was selected", fill=FALSE)
        if(verbose) cat("with threshold knot", round(zetas[mystep],3), " !", fill=TRUE)
        
        ## Terminate if all terminal nodes are length 2 or smaller.
        too.short = unlist(lapply(Tt[[mystep]], function(mypair){Ecurr[mypair[1],mypair[2]] - Scurr[mypair[1], mypair[2]] <=1}))
        if(all(too.short)){
            if(verbose) cat("Ended early, at step", mystep, fill=TRUE)
            break;
        }
    }
    ## END of Main loop
    G = trim(G,"row")
    
    ## trim and return things
    return(list(S = trim(S),
                E = trim(E),
                A = trim(A),
                T = trim(Tt),
                B = trim(B),
                Z = trim(Z),
                G = G,
                u = rep(0,nrow(G)),
                zetas = trim(zetas),
                numsteps = mystep
                ))
}

