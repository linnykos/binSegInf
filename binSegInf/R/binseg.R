##' Main function for binary segmentation for fixed threshold. This is actually
##' a wrapper for binary segmentation with fixed threshold. It creates an
##' environment and creates the variables there, then runs
##' binseg.by.thresh.inner() all in this environment, and returns the relevant
##' guy
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param j The depth of the recursion on hand.
##' @param k The indexing of the node location, from left to right, in the
##'     \emph{complete} binary tree.
##' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
##'     for the recursion.
##' @param y The original data.
##' @param verbose Set to true if you'd like to see algorithm progression.
##' @param return.env Set to true if you'd like to get the environment
##'     containing the algorithm output, instead of a list.
##' @import pryr

binseg.by.thresh = function(y, thresh, s=1, e=length(y), j=0, k=1, verbose=FALSE, return.env=FALSE){
    n = length(y)
    if(n > 20) 
{
        stop(paste0("You'll use up a /lot/ of memory, so I'm stopping you.",
                    "This is my fault, not yours. I'm using giant mostly empty matrices. ",
                    "I'm embarassed. I'll change this soon. Prod me if you need me to do it now."))
        }

    ## Create new environment |env|
    env = new.env()
    env$slist = env$elist = env$blist = env$Blist = env$zlist = env$Zlist =
        ## matrix(NA, nrow = n, ncol = 2 ^(n/2))
        ## matrix(NA,nrow=n,ncol=3)
        cplist(n)
    env$y = y

    ## Run binary segmentation on |env|
    binseg.by.thresh.inner(y, thresh, s, e, j, k, verbose, env=env)

    ## Gather output from |env| and return it.
    bs.output = list(slist = trim(env$slist),
                     elist = trim(env$elist),
                     blist = trim(env$blist),
                     Blist = trim(env$Blist),
                     zlist = trim(env$zlist),
                     Zlist = trim(env$Zlist),
                     cp = trim(env$blist$mat[,3]),
                     cp.sign = trim(env$zlist$mat[,3]),
                     y = y,
                     thresh = thresh)
    if(return.env){
        return(env)
    } else{
        return(bs.output)
    }
}





##' Inner function for binary segmentation with fixed threshold. The wrapper
##' \code{binseg.by.thresh()} is intended to be used by user. Note, when
##' \code{thresh} is set to zero, then this can be used to collect the
##' unbalanced haar wavelet basis.
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param j The depth of the recursion on hand.
##' @param k The indexing of the node location, from left to right, in the
##'     \emph{complete} binary tree.
##' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
##'     for the recursion.
##' @param y The original data.
##' @param n The length of the data \code{y}.

binseg.by.thresh.inner <- function(y, thresh, s=1, e=length(y), j=0, k=1, verbose=F, env=NULL){

    n = length(y)
    
    ## If segment is 2-lengthed, terminate
    if(e-s<=1){
        env$slist = add(env$slist,j,k,s)
        env$elist = add(env$elist,j,k,e)
        return() 
    ## Otherwise, calculate CUSUMs
    } else {
        all.bs = (s:(e-1))
        all.cusums = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
        names(all.cusums) = all.bs
        all.cusums = c(rep(NA,s-1),all.cusums, rep(NA,n-s))
        sn = sign(all.cusums) 
        
        ## Obtain breakpoint and its sign
        b = which.max(abs(all.cusums))
        z = sn[b]
        if(verbose) cat("the changepoint", b, "is selected", fill=T)
        
        ## Check threshold exceedance, then store
        if(abs(all.cusums[b]) < thresh){
            env$Blist = add(env$Blist,j+1,k,b)
            env$Zlist = add(env$Zlist,j+1,k,z)
            env$slist = add(env$slist,j,  k,s)
            env$elist = add(env$elist,j,  k,e)
            return(env)
        } else { 
            if(verbose) cat("the biggest cusum was", all.cusums[b], "which passed the threshold:", thresh,  fill=T)

            env$blist = add(env$blist, j+1,k, b)
            env$Blist = add(env$Blist, j+1,k, b)
            env$zlist = add(env$zlist, j+1,k, z)
            env$Zlist = add(env$Zlist, j+1,k, z)
            env$slist = add(env$slist, j  ,k, s)
            env$elist = add(env$elist, j  ,k, e)
        }
        
        ## Recurse
        binseg.by.thresh.inner(y, thresh, s, b, j+1, 2*k-1, verbose, env=env)
        binseg.by.thresh.inner(y, thresh, b+1, e, j+1, 2*k, verbose, env=env)
    }
}


##' Function to carry out fixed-number-of-steps binary segmentation.
##' @param y numeric vector, data
##' @param numsteps desired number of changepoints
##' @param verbose set to true for algorithm run details. 
##' @import Matrix
##' @export

binseg.by.size = function(y,numsteps,verbose=FALSE){

    ## Basic checks
    if(numsteps > length(y)-1) stop(paste("You should ask for less than", length(y), "steps!"))
    
    ## Initialize things
    B = Z = rep(NA,length(y)) 
    Bcurr = Ccurr = Zcurr = Scurr = Ecurr = cplist(numsteps+1)

    S = E = A = Tt = Tcurr = Acurr = lapply(1:length(y),function(i) c())
    Scurr = add(Scurr,1,1,1)
    Ecurr = add(Ecurr,1,1,length(y))
    Tcurr[[1]] = c(1,1)
    Acurr[[1]] = c()
    
    jk  = list()
    G = matrix(NA, ncol = length(y), nrow = 2*length(y)*numsteps)
    G.jk = cplist(numsteps+1) ## G.jk will store the row numbers corresponding
                              ## to the sign characterizers of (b,j,k) over all
                              ## b, for each (j,k)
    Gn = 0 

    Tcurr.which.new = 1
    
    ## Main loop
    for(mystep in 1:numsteps){
        if(verbose) cat("At step", mystep, " ", fill=TRUE)

        ## Get all candidate changepoints
        Gn.beginning.of.step = Gn

        ## Go through non-null terminal nodes, but only calculate new things in Tcurr!
        for(ii in Tcurr.which.new){
            j = Tcurr[[ii]][1]
            k = Tcurr[[ii]][2]

            cusums = getcusums(s = extract(Scurr,j,k),
                               e = extract(Ecurr,j,k),
                               y = y)
            ## Characterize signs
            signed.cusummat = (cusums$contrasts) * (cusums$signs)
            G[(Gn+1):(Gn+nrow(signed.cusummat)),] = signed.cusummat 
            G.jk = add(G.jk, j, k, paste0(Gn+1,":",Gn+nrow(signed.cusummat)))
            Gn = Gn+nrow(signed.cusummat)

            ## Find cusum maximizer
            Bcurr = add(Bcurr,j,k,cusums$bmax)
            Ccurr = add(Ccurr,j,k,cusums$cusum)
            Zcurr = add(Zcurr,j,k,cusums$signs[cusums$bmax.cusums])
        }

        ## Get all maximizers of each block (Note: one /block/ of indices per node)
        all.rownums.this.step = unlist(lapply(rid.null(Tcurr),function(jk){mydeparse(extract(G.jk, jk[1], jk[2]))}))
        all.cusums.this.step = lapply(rid.null(Tcurr),function(jk){extract(Ccurr, jk[1], jk[2])})
        all.cusums.this.step = unlist(all.cusums.this.step)
        max.jk.this.step = rid.null(Tcurr)[[which.max(all.cusums.this.step)]]

        ## Get the block that contains the maximizer, and find the maximum from that block
        max.block.rownums = mydeparse(extract(G.jk, max.jk.this.step[1], max.jk.this.step[2]))
        max.block.cusums = G[max.block.rownums,]%*% cbind(y)
        max.rownum = max.block.rownums[which.max(max.block.cusums)]
        
        ## Characterize cusum-maximizing rows /of this step/
        all.rownums.this.step = all.rownums.this.step[all.rownums.this.step!=max.rownum]  
        all.rownums.this.step = unlist(lapply(rid.null(Tcurr),function(jk){mydeparse(extract(G.jk, jk[1], jk[2]))}))
        this.step.G = G[all.rownums.this.step[all.rownums.this.step != max.rownum],,drop=FALSE]
        comparison.cusummat = t(apply(this.step.G, 1, function(myrow){ G[max.rownum,,drop=FALSE] - myrow}))

        ## Add it to G
        G[(Gn+1):(Gn+nrow(comparison.cusummat)),] = comparison.cusummat
        Gn = Gn+nrow(comparison.cusummat)

        ## Record knot as CUSUM maximizer
        jmax = max.jk.this.step[1]
        kmax = max.jk.this.step[2]
        
        ## Update terminal and active node set
        which.duplicate = which(sapply(Tcurr, function(myjk){all.equal(myjk, c(jmax, kmax))==TRUE}))
        Tcurr[[which.duplicate]] <- c(jmax + 1, 2*kmax - 1)
        Tcurr[[mystep+1]] <- c(jmax + 1, 2*kmax)
        Tcurr.which.new = c(which.duplicate, mystep+1)
        Acurr[[mystep+1]] <- c(jmax,kmax)
        
        ## Update Scurr and Ecurr for the /new/ nodes
        ## if(mystep==6) browser()
        Scurr = add(Scurr, jmax+1, 2*kmax-1, extract(Scurr,jmax,kmax))
        Ecurr = add(Ecurr, jmax+1, 2*kmax-1, extract(Bcurr,jmax,kmax))
        Scurr = add(Scurr, jmax+1, 2*kmax, extract(Bcurr,jmax,kmax) + 1)
        Ecurr = add(Ecurr, jmax+1, 2*kmax, extract(Ecurr,jmax,kmax))

        ## Rid Tcurr of all nodes that represent 2-length segments
        riddance.of.short.segments = rid_jk_nicely_from_Tcurr(Tcurr, Scurr, Ecurr, Tcurr.which.new)
        Tcurr = riddance.of.short.segments$Tcurr
        Tcurr.which.new = riddance.of.short.segments$Tcurr.which.new
        
        ## Take snapshot
        S[[mystep]] = trim(Scurr)
        E[[mystep]] = trim(Ecurr)
        A[[mystep]] = trim(Acurr)
        Tt[[mystep]] = rid.null(Tcurr)
        B[mystep] = extract(Bcurr, jmax, kmax)
        Z[mystep] = extract(Zcurr, jmax, kmax)
        jk[[mystep]] = c(jmax,kmax) 
        
        if(verbose) cat("From candidates", extract(Scurr,jmax,kmax),  ":",extract(Ecurr, jmax, kmax), " ")
        if(verbose) cat("breakpoint", extract(Bcurr,jmax,kmax), "was selected", fill=TRUE)
        
        ## Terminate if all terminal nodes are length 2 or smaller.
        too.short = unlist(lapply(Tt[[mystep]], function(mypair){extract(Ecurr, mypair[1],mypair[2]) -
                                                                     extract(Scurr, mypair[1], mypair[2]) <=1}))
        if(all(too.short)){
            if(verbose) cat("Ended early, at step", mystep, fill=TRUE)
            break;
        }
    }
    ## END of Main loop
    G = trim(G,"row")

    ## Trim and return things
    return(list(S = trim(S),
                E = trim(E),
                A = trim(A),
                T = trim(Tt),
                B = trim(B),
                Z = trim(Z),
                G = G,
                u = rep(0,nrow(G)),
                numsteps = mystep
                ))
}
