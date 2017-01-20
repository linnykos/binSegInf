##genlassoinf can be found at https://github.com/robohyun66/genlassoinf/tree/master/code

context("Test flasso implementation against genlassoinf")

dualpathSvd2 <- function(y, D, approx=FALSE, maxsteps=2000, minlam=0,
                         rtol=1e-7, btol=1e-7, verbose=FALSE, object=NULL,
                         ctol=1e-10, cdtol=1e-4, do.gc=F){
  
  # Error checking
  stopifnot(ncol(D) == length(y))
  
  nk = 0
  ss = list() # list of ss
  tab = matrix(NA,nrow=maxsteps,ncol=6) # list of where
  colnames(tab) = rev(c("after hit vs leave (leave wins)",
                        "after hit vs leave (hit wins)",
                        "after leave times",
                        "after leave eligibility (c<0)",
                        "after hitting event",
                        "after viable hits"))
  
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)
    
    # Initialize Gamma matrix (work in progress)
    # Gammat = Matrix(NA,nrow= maxsteps*ncol(D)*4 ,ncol=n)
    
    # Compute the dual solution at infinity, and
    # find the first critical point
    In = diag(1,n)
    sv = svdsolve(t(D),y,rtol)
    uhat = as.numeric(sv$x[,1])        # Dual solution
    q = sv$q                           # Rank of D
    
    ihit = which.max(abs(uhat))   # Hitting coordinate
    hit = abs(uhat[ihit])         # Critical lambda
    s = Sign(uhat[ihit])          # Sign
    k = 1
    ss[[k]] = s 
    
    if (verbose) {
      cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                  hit,ihit,1))
    }
    
    # Now iteratively find the new dual solution, and
    # the next critical point
    
    # Things to keep track of, and return at the end
    buf = min(maxsteps,1000)
    u = matrix(0,m,buf)        # Dual solutions
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom
    action = numeric(buf)      # Action taken
    upol = c()                 # Constant in polyhedral constraint
    
    lams[1] = hit
    action[1] = ihit
    h[1] = TRUE
    df[1] = n-q
    u[,1] = uhat
    
    
    # add rows to Gamma
    tDinv = MASS::ginv(as.matrix(t(D)))
    
    # rows to add, for first hitting time (no need for sign-eligibility--just add all sign pairs)
    #Gammat = diag(Sign(uhat)  %*% tDinv)  
    M = matrix(s*tDinv[ihit,], nrow(tDinv[-ihit,]), n, byrow=TRUE) 
    Gammat = rbind(M + tDinv[-ihit,], 
                   M - tDinv[-ihit,])
    
    tab[k,2] = nrow(Gammat)
    
    nk = nrow(Gammat)
    
    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    Ds = D[ihit,]*s            # Vector t(D[B,])%*%s
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
    
    
  } else {
    # If iterating an already started path
    # Grab variables needed to construct the path
    lambda = NULL
    for (j in 1:length(object)) {
      if (names(object)[j] != "pathobjs") {
        assign(names(object)[j], object[[j]])
      }
    }
    for (j in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[j], object$pathobjs[[j]])
    }
    lams = lambda
  }
  
  tryCatch({
    while (k<=maxsteps && lams[k-1]>=minlam) {
      
      if(do.gc) gc()
      if(verbose){ cat('\n'); show.glutton(environment(),4)}
      
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        action = c(action,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
      }
      
      ##########
      # If the interior is empty, then nothing will hit
      if (r==m) {
        a = b = numeric(0)
        hit = 0
        q = 0
      }
      # Otherwise, find the next hitting time
      else {
        In = diag(1,n)
        sv = svdsolve(t(D1),cbind(y,Ds),rtol) #sv = svdsolve(t(D1),cbind(y,Ds,In),rtol)
        a = as.numeric(sv$x[,1])  # formerly a = as.numeric(D3 %*% y)
        b = as.numeric(sv$x[,2])
        D3 = MASS::ginv(as.matrix(t(D1)))# formerly as.matrix(sv$x[,3:(n+2)]) 
        # meant to be pseudoinverse of t(D[-I,])
        
        q = sv$q
        shits = Sign(a)
        hits = a/(b+shits);
        
        
        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]+btol] = 0
        hits[hits>lams[k-1]] = lams[k-1]
        
        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]
        
        # Gamma Matrix!
        # rows to add, for viable hitting signs:
        tDinv = D3
        rows.to.add = (if(length(shits)>1){
          do.call(rbind, lapply(1:length(shits), function(ii){shits[ii] * tDinv[ii,]  }))
        } else {
          shits* tDinv
        })
        
        Gammat = rbind(Gammat, rows.to.add)
        tab[k,1] = nrow(Gammat)
        
        # rows to add, for hitting event: (This is just dividing each row of tDinv by corresponding element of tDinv%*%Ds+shit)
        A = asrowmat(D3/(b+shits)) # tDinv / as.numeric(tDinv %*% Ds + shits)
        if(nrow(A)!=1){
          nleft = nrow(A[-ihit,])
          if(is.null(nleft)) nleft = 1
          M = matrix(A[ihit,], nrow = nleft, ncol = n, byrow = TRUE) 
          Gammat = rbind(Gammat, M - A[-ihit,])
        }
        tab[k,2] = nrow(Gammat)
      }
      
      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }
      
      # Otherwise, find the next leaving time
      else {
        c = as.matrix(s*(D2%*%(y-t(D1)%*%a)))
        d = as.matrix(s*(D2%*%(Ds-t(D1)%*%b)))
        
        
        # round small values of c to zero (but not d)
        #cdtol = 1E-10
        c[abs(c) < cdtol] = 0
        
        # get leave times
        leaves = c/d
        
        # identify which on boundary set are c<0 and d<0
        Ci = (c < 0)
        Di = (d < 0)
        Ci[!B] = Di[!B] = FALSE
        CDi = (Ci & Di)
        
        # c and d must be negative at all coordinates to be considered
        leaves[c>=0|d>=0] = 0
        
        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        super.lambda = leaves>lams[k-1]+btol
        closeto.lambda = (lams[k-1] < leaves ) & (leaves  < lams[k-1]+btol)
        leaves[leaves>lams[k-1]+btol] = 0
        leaves[leaves>lams[k-1]] = lams[k-1]
        
        # If a variable just entered, then make sure it 
        # cannot leave (added from lasso.R)
        if (action[k-1]>0) leaves[r] = 0
        
        # index of leaving coordinate
        ileave = which.max(leaves)
        leave = leaves[ileave]
        
        # Gamma Matrix!!
        # rows to add, for leaving event:
        if(dim(D1)[1]==0) D1 = rep(0,ncol(D1)) # temporarily added because of 
        # dimension problem in next line, 
        # at last step of algorithm
        gmat = s*(D2%*%(In - t(D1)%*%D3)) # coefficient matrix to c
        
        # close-to-zero replacement is hard-coded in
        gmat[abs(c)<cdtol,] = rep(0,ncol(gmat))
        
        # we still want to see that gmat&%y ~= c
        if(!(max(gmat%*%y-c) < ctol)) print(max(gmat%*%y-c))
        
        gd = gmat / as.numeric(d)
        
        
        # hard coding in the zero replacements, to make it identical with lea
        gd[c>=0,] = rep(0, ncol(gd))        # re-doing zero-replacement in the g/d matrix
        gd[super.lambda,] = rep(0,ncol(gd)) # re-doing larger-than-lambda-replacement 
        #gd[closeto.lambda,] = rep(0,ncol(gd)) # re-doing close-to-lambda-replacement (not sure how to make the i'th leaving time gd[i,]%*%y == lam[k-1] properly; solving to get some gd[i,] is like finding a non-unique solution to an overdetermined system; because such a gd[i,] is not unique, how do I know that adding this row to Gamma won't do weird and mysterious things?)
        
        #if( (length(Di)!=0) & (which(closeto.lambda) %in% which(Di))) print("closeto.lambda replacement SHOULD have happenned (but didn't).")
        
        # add rows that ensure c<0 #(only in )
        Gammat <- rbind(Gammat, gmat[Ci&Di,]*(-1),
                        gmat[(!Ci)&Di,])
        #        print("after leave eligibility (c<0)")
        #        print(nrow(Gammat))
        tab[k,3] = nrow(Gammat)
        
        # get rid of NA rows in Gammat (temporary fix)
        missing.rows = apply(Gammat, 1, function(row) any(is.na(row)))
        if(sum(missing.rows)>=1){ Gammat <- Gammat[-which(missing.rows),] }
        
        # add rows for maximizer
        CDi = (Ci & Di)
        CDi[ileave] = FALSE
        CDind = which(CDi)
        
        Gammat <- rbind(Gammat, gd[rep(ileave,length(CDind)),] - gd[CDind,])
        #        print("after leave times")
        #        print(nrow(Gammat))
        tab[k,4] = nrow(Gammat)
      }
      ##########
      # Stop if the next critical point is negative
      
      if (hit<=0 && leave<=0) {break}
      
      # If a hitting time comes next
      if (hit > leave) {
        
        # Record the critical lambda and solution
        lams[k] = hit
        action[k] = I[ihit]
        h[k] = TRUE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        u[,k] = uhat
        
        
        # add row to Gamma to characterize the hit coming next
        if(!approx)  {
          # this is literally h_k - l_k > 0
          Gammat = rbind(Gammat,  A[ihit,] - gd[ileave,])
          #          print("after hit vs leave (hit wins)")
          #          print(nrow(Gammat))
          tab[k,5] = nrow(Gammat)
        }
        
        nk = c(nk,nrow(Gammat))
        
        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        Ds = Ds + D1[ihit,]*shit
        s = c(s,shit)
        D2 = rbind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
        ss[[k]] = s
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }
      
      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and solution
        lams[k] = leave
        action[k] = -B[ileave]
        h[k] = FALSE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        u[,k] = uhat
        
        
        # add row to Gamma to characterize the leave coming next
        if(!approx)  {
          Gammat = rbind(Gammat, - A[ihit,] + gd[ileave,])
          tab[k,5] = nrow(Gammat)
        }
        nk = c(nk,nrow(Gammat))
        
        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        Ds = Ds - D2[ileave,]*s[ileave]
        s = s[-ileave]
        D1 = rbind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]
        ss[[k]] = s
        
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
        
      }
      
      # Step counter + other stuff
      k = k+1
      # resetting tDinv
      #tDinv = t(as.matrix(rep(NA,length(tDinv)))[,-1,drop=F])
    }
  }, error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
                        " partial path is being returned.)",sep="")
    warning(err)})
  
  # Trim
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  u = u[,Seq(1,k-1),drop=FALSE]
  
  # Save needed elements for continuing the path
  pathobjs = list(type="svd", r=r, B=B, I=I, approx=approx,
                  k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, h=h,
                  rtol=rtol, btol=btol, s=s, y=y)
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  } else if (lams[k-1]<minlam) {
    # If we reached the minimum lambda
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  } else {
    # Otherwise, note that we completed the path
    completepath = TRUE
  }
  if (verbose) cat("\n")
  
  colnames(u) = as.character(round(lams,3))
  
  beta <-  apply(t(D)%*%u,2,function(column){y-column}) 
  ss = c(NA,ss)
  states = get.states(action)
  
  return(list(lambda=lams,beta=beta,fit=beta,u=u,hit=h,df=df,y=y,ss=ss,states=states,
              completepath=completepath,bls=y,pathobjs=pathobjs, 
              Gammat=Gammat, nk = nk, action=action,tab=tab))
}


# general function that makes various D matrices
makeDmat = function(m, type = c("trend.filtering","2Dgraph"), order=0){
  
  type = match.arg(type)
  D = NA
  if(type == "trend.filtering"){ 
    # handles fused lasso and more
    D = dual1d_Dmat(m)
    if(order>=1){
      for(jj in 1:order){ D = dual1d_Dmat(m-jj) %*% D }
    }
    return(D)  
  } else {
    stop("Not coded yet!")
  }
}

dual1d_Dmat = function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}

getGammat.naive = function(obj, y, condition.step){#, usage = c("fl1d", "dualpathSvd")){
  
  n = length(y)
  
  # Error checking
  #  usage <- match.arg(usage)
  if(length(obj$action) < condition.step ) stop("\n You must ask for polyhedron from less than ", length(obj$action), " steps! \n")
  #  if(length(obj) < 11 && "completepath" %in% objects(obj)) stop("|obj| and |usage| are mismatched")
  
  # get appropriate sized Gamma matrix from selection step |condition.step|
  #  if(usage == "fl1d"){
  #    myGammat = obj$Gammat[1:obj$nk[condition.step],]
  #  } else {
  myGammat = obj$Gamma[1:obj$nk[condition.step],]
  return(list(Gammat = myGammat,condition.step=condition.step))
}




#' get k'th dik vector, for the segment or spike test.
getdvec = function(obj, y=NULL, k, klater = k, type =c("spike","segment","hybridsegment","general"), n){
  # obj : either output from fl1d() or from dualpathSvd2()
  # y     : data vector (same as the one used to produce path)
  # k     : algorithm's step to use.
  # klater: later step to condition on
  
  type <- match.arg(type)
  if(k > klater) stop("Attempting to form contrast from a later step than the conditioning step.")
  
  
  # ik and sk are the index and sign of the jump selected at step k
  
  ik = (obj$pathobjs)$B[k]
  sk = (obj$pathobjs)$s[k] 
  breaks = (obj$pathobjs)$B
  
  if(type == "spike"){
    
    dik = rep(0,length(y))
    dik[ik] = -1
    dik[ik+1] = 1
    dik = sk * dik     # this ensures that the gap we're testing is positive. i.e. the jump sk(y_(ik) - y_(ik-1)) always positive!
    
  } else if (type == "segment"){
    
    ## extract usual segment test endpoints
    Ks = makesegment(breaks=breaks,k=k,klater=klater,n=length(y))
    K = Ks$K
    Kmin = Ks$Kmin
    Kmax = Ks$Kmax
    
    # form vector
    dik = rep(0,length(y))    
    dik[Kmin:K] <- (Kmax - K)/(Kmax - Kmin + 1)
    dik[(K+1):Kmax] <- -(K - Kmin + 1)/(Kmax - Kmin + 1)
    dik <- -sk *dik
    
  } else if (type == "hybridsegment"){
    ## make sure you run it to a reasonable (or full set of) steps
    
    ## extract minimum length
    sortedbreaks = sort(breaks)
    if(length(sortedbreaks)<=1){
      lens = 0
    } else {
      lens = sortedbreaks [2:length(sortedbreaks )] - sortedbreaks [1:(length(sortedbreaks )-1)]
    }
    minlen = mean(lens)*1.5 
    
    ## extract usual segment test endpoints
    Ks = makesegment(breaks=breaks,k=k,klater=klater,n=length(y))
    K = Ks$K
    Kmin = Ks$Kmin
    Kmax = Ks$Kmax
    
    ## replace Kmax and Kmin if necessary
    if( (K - Kmin) > minlen ) Kmin = ceiling(K - minlen)
    if( (Kmax - K) > minlen ) Kmax = floor(K + minlen)
    
    # form vector
    dik = rep(0,length(y))
    dik[Kmin:K] <- (Kmax - K)/(Kmax - Kmin + 1)
    dik[(K+1):Kmax] <- -(K - Kmin + 1)/(Kmax - Kmin + 1)
    dik <- -sk *dik
    
  } else if (type == "general") {     
    stop("Not coded yet!")  
  } else {
    stop("Not coded yet!")
  }
  return(dik)
}


makesegment  = function(breaks, k, klater, n){
  
  if(length(breaks)<k) stop("not enough breaks!! k > number of breaks")
  K <- breaks[k] # is the index of the jump selected at step k
  
  # whether or not to condition on a later step (|klater|)
  kk <- klater
  
  relevantbreaks = (if(kk==1) c() else breaks[1:kk])
  endpoints = c(1,n)
  allbreaks <- c(endpoints, relevantbreaks)
  allbreaks <- sort(allbreaks)
  
  if(K %in% allbreaks) allbreaks = allbreaks[-which(allbreaks == K)]
  allbreaks = sort(unique(c(allbreaks, endpoints))) #temporary fix just in case the global endpoints are detected..
  min.index <- max(sum(allbreaks< K),1)             #temporary fix continued 
  
  Kmin <- allbreaks[min.index]
  Kmax <- allbreaks[min.index + 1]    
  
  if(Kmin != 1) Kmin = Kmin + 1 # special case handling
  
  return(list(Kmin=Kmin,K=K,Kmax=Kmax))
}

pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
}


poly.pval <- function(y, G, u, v, sigma, bits=NULL) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  pv = tnorm.surv(z,0,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

# Main confidence interval function

poly.int <- function(y, G, u, v, sigma, alpha, gridrange=c(-100,100),
                     gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {
  
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }
  
  int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  tailarea = c(fun(int[1]),1-fun(int[2]))
  
  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea))
}


tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 
  
  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

##' Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, using
##' multi precision floating point calculations thanks to the Rmpfr package
mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  
  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

svdsolve <- function(A,b, rtol=1e-7) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii),s=s))
}

Sign <- function(x) {
  return(-1+2*(x>=0))
}

Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

get.states = function(action.obj){
  
  # Obtain a list of actions at each step.
  actionlists = sapply(1:length(action.obj),function(ii) action.obj[1:ii] )
  
  # Helper function to extract the final state after running through (a list of) actions
  get.final.state = function(actionlist){
    if(length(actionlist)==1) return(actionlist)
    to.delete = c()
    for(abs.coord in unique(abs(actionlist))){
      all.coord.actions = which(abs(actionlist)==abs.coord)
      if( length(all.coord.actions) > 1 ){
        if( length(all.coord.actions) %%2 ==1){
          to.delete = c(to.delete, all.coord.actions[1:(length(all.coord.actions)-1)])
        } else {
          to.delete = c(to.delete, all.coord.actions)
        }
      }
    }
    if(!is.null(to.delete)){
      return(actionlist[-to.delete])
    } else {
      return(actionlist)
    }
  }
  states = lapply(actionlists, get.final.state)
  states = c(NA,states)
  
  return(states)
}

asrowmat = function(obj){
  objmat = as.matrix(obj)
  if(ncol(objmat)==1 & ncol(objmat) < nrow(objmat)){
    objmat = t(objmat)
  }
  return(objmat)
}


################################################


test_that("dimension of gamma is the same", {
  set.seed(10)
  #use justin's code
  n <- 10
  y <- rnorm(n, 0, 1)
  D <- makeDmat(n, ord = 0)
  mypath <- dualpathSvd2(y, D, 1, approx = T)
  G <- getGammat.naive(obj = mypath, y = y, condition.step = 1)
  d <- getdvec(mypath, y, 1, 1, type = "segment")
  pval <- pval.fl1d(y, G$Gammat, d, 1)
  
  #use our code
  obj <- fLasso_fixedSteps(y, 1)
  res <- polyhedra(obj)
  
  expect_true(all(dim(G$Gammat) == dim(res$gamma)))
})

test_that("dimension of gamma is the same", {
  set.seed(10)
  #use justin's code
  n <- 10
  y <- rnorm(n, 0, 1)
  D <- makeDmat(n, ord = 0)
  mypath <- dualpathSvd2(y, D, 1, approx = T)
  G <- getGammat.naive(obj = mypath, y = y, condition.step = 1)
  d <- getdvec(mypath, y, 1, 1, type = "segment")
  pval <- pval.fl1d(y, G$Gammat, d, 1) 
  
  #use our code
  obj <- fLasso_fixedSteps(y, 1)
  res <- polyhedra(obj)
  
  match_vec <- numeric(nrow(res$gamma))
  match_value <- numeric(nrow(res$gamma))
  for(i in 1:nrow(res$gamma)){
    vec <- sapply(1:nrow(G$Gammat), function(x){
      sqrt(sum((G$Gammat[x,] - res$gamma[i,])^2))
    })
    
    match_vec[i] <- which.min(vec)
    match_value[i] <- min(vec)
  }
  
  expect_true(length(unique(match_vec)) == length(match_vec))
  expect_true(all(match_value < 1e-4))
})


