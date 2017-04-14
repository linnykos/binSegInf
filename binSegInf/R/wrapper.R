##' Embed stopping time into the path.
##' @param obj A list object containing an element named \code{"cp"}, which is a
##'     integer vector with changepoints, in the order they were selected by the
##'     sequential algorithm
##' @param sigma Standard deviation of noise generating the data, in \code{y}
##' @param type One of \code{c("bic","ebic","aic")}. Defaults to
##'     \code{"bic"}. Only BIC is possible for now.
##' @param consec Number of rises in IC before stop is to happen. Defaults to 2.
##' @return list containing two objects: ''poly'' and ''stoptime''.
##' @example examples/ic_wrapper-example.R
##' 
##' @export
ic_wrapper <- function(obj, y, consec=2, maxsteps=length(cp), sigma,
                       type = "bic"){## Basic checks
    stopifnot(c("cp") %in% objects(obj))
    if(type!="bic") stop("Only BIC is possible, for now.")
    
    ## Get ic information
    tryCatch({
        ic_obj = get_ic(obj$cp, y, consec=2, sigma=sigma, type = type)
    },
    warning = function() return(NULL)
    ) 
    newpoly = ic_to_poly(ic_obj)
    
    ## Return 
    return(list(poly=newpoly, stoptime = ic_obj$stoptime+1))
}

##' Takes IC vector and creates polyhedra.
##' @param resid list of residual vectors
##' @param obj object of class \code{ic} from \code{get_ic()}, which contains
##'     things needed for sequential information criteria (ic) comparison.
##' 
##' @return Object of class \code{polyhedra}, for sequential IC comparisons.
ic_to_poly <- function(obj){

    ## Basic checks
    stopifnot(is_valid.ic(obj))

    ## Get order of ICs
    seqdirs = c(.getorder(obj$ic))

    ## Make empty things
    newrows = matrix(NA, nrow = 2*(obj$stoptime+obj$consec),
                     ncol = length(obj$y))
    newu = rep(NA, 2*(obj$stoptime+obj$consec))
    irow = 0
    
    ## Collect halfspaces
    for(jj in 1:(obj$stoptime + obj$consec)){

    
        residual = obj$resid[[jj+1]]
        const    = obj$pen[jj+1] - obj$pen[jj]

        if(seqdirs[jj+1] > 0){
            ## Add one row \sqrt{C} < z_a \times a^Ty 
            newrows[irow+1,] = sign(t(residual)%*%obj$y) * residual/sqrt(sum(residual^2))
            newu[irow+1] = sqrt(const)
            irow = irow + 1
        } else {
            ## Add two rows -\sqrt{C} < z_a \times a^Ty < \sqrt{C}
            newrows[irow+1,] = sign(t(residual)%*%obj$y) * residual/sqrt(sum(residual^2))
            newrows[irow+2,] = -sign(t(residual)%*%obj$y) * residual/sqrt(sum(residual^2))
            newu[irow+1] = -sqrt(const) 
            newu[irow+2] = -sqrt(const)
            irow = irow + 2
        }
    }
    ## Form polyhedra
    poly = polyhedra(obj = trim(newrows), u = trim(newu))
    return(poly)
}


##' Function that takes changepoint location vector \code{cp} (in the order that
##' they entered, from fixed step SBS of fixed step WBS), and calcualtes BIC of
##' each intermediate changepoint model (piecewise constant underlying mean
##' Gaussian model).
##' @param cp Vector of changepoints
##' @param y Data vector.
##' @param sigma Noise standard deviation.
##' @param verbose \code{TRUE} to print progress
##'
##' @return A list containing two elements: a numeric vector of BICs, and a list
##'     containing basis vectors \eqn{a} of the projection \eqn{P=aa^T}from one
##'     changepoint column basis space to the next.
get_ic <- function(cp, y, sigma, consec=2, maxsteps=length(cp), type="bic", verbose=FALSE){
    ## Basic checks
    if(type!="bic") stop("Only BIC is coded so far!")
    
    ## Collect things
    n = length(y)
    D = dual1d_Dmat(n)
    ic = pen = RSS = rep(NA, maxsteps)
    resid = list()
    
    ## Collect BIC at each step 0 ~ (maxsteps-1)
    for(ii in 1:pmin(maxsteps,length(cp))){
        ## if(verbose)  cat('step', ii, '\n')
         cat('step', ii, '\n')
        
        ## Form proj null(D_{-B}) by orth proj onto row(D_{-B}) = col(t(D_{-B})) ~= tD
        tD = cbind(t(D)[,-cp[1:ii]])
        rr = rankMatrix(tD)
        tDb = svd(tD)$u[,1:rr]
        curr.proj = .proj(tDb)
        y.fitted = (diag(1,n) - curr.proj) %*% y
        
        ## Obtain RSS and penalty
        myRSS = sum( (y - y.fitted)^2 )
        mydf  = n-rr
        if(ii==1) prev.df = mydf
        mypen = (sigma^2) * mydf * log(n) 
        
        ## Obtain (2norm-scaled) residual projection vectors
        if(ii==1){
            myresid = rep(NA,n)
        } else {
            myresid = svd(curr.proj - prev.proj)$u[,1]
            myresid = myresid / sqrt(sum((myresid)^2))
        }
        
        ## Store BIC and resid proj vector
        ic[ii] <- myRSS + mypen
        pen[ii] <- mypen
        RSS[ii] <- myRSS
        resid[[ii]] <- myresid
        
        ## Update things for next step.
        prev.proj = curr.proj
    }

    ## Obtain stoptime
    stoptime = .whichrise(ic,consec) - 1 
    stoptime = pmin(stoptime, length(y)-consec-1)

    ## Return NULL if path hasn't stopped
    if(!(stoptime+consec < maxsteps)){
        warning(paste('IC rule using', consec, 'rises hasnt stopped!'))
    }
    if(stoptime==0) warning('Stoptime is zero!')

    obj = structure(list(ic=ic, consec=consec, resid=resid, pen=pen, RSS=RSS,
                         stoptime=stoptime,y=y, type=type), class="ic")
    return(obj)
}


##' is_valid for ic objects
##' @return TRUE if valid
is_valid.ic <- function(obj){
    if(!(all(names(obj)%in% c("ic","consec","resid","pen","type","stoptime","RSS",
                              "y") ))){
        stop("obj needs to be produced from get_ic()")
    }
    TRUE
}



##' Helper function to project on /column/ space of matrix.
.proj <- function(mymat){
    return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))
}

##' Returns a sequence of +1 and -1 for sequential incr and decrements in a
##' vector; assume first step always dips; _almost_ always true
##' @param ic a numeric vector
.getorder = function(ic){
  return(c(NA,sign(ic[2:(length(ic))] - ic[1:(length(ic)-1)])))
}




##' Locates first point of \code{consec} rises in IC path
##' @param ic numeric vector containing information criteria.
##'
##' @return First point at which ic rises.
##' @export
.whichrise = function(ic, consec = 2, direction=c("forward","backward")){

  direction = match.arg(direction)
  if(direction != "forward") stop("That direction IC selection is not coded yet.")
  if(length(ic) < consec+1) stop("Not enough steps to do forward sequential BIC/AIC!")  

  ind = 1
  done = FALSE
  while(ind < (length(ic)-consec+1) ){
    ictol = 1E-10
    if( all(ic[(ind+1):(ind+consec)] > ic[(ind):(ind+consec-1)] + ictol )) break
    ind = ind+1
  }
  return(pmin(pmax(ind,1),length(ic)))
}
