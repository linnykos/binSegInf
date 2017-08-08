#' Generate polyhedra matrix from bsFs object
#'
#' Forms both Gamma matrix and u vector
#'
#' @param obj bsFs object
#' @param ... not used
#'
#' @return An object of class polyhedra
#' @export
polyhedra.bsFs <- function(obj, ...){
  is_valid(obj)

  n <- .get_startEnd(obj$tree$name)[2]
  numSteps <- obj$numSteps
  comp.lis <- .list_comparison(obj)
  sign.vec <- sign(jump_cusum(obj))
  gamma.row.lis <- vector("list", numSteps)

  for(i in 1:numSteps){
    losing.mat <- comp.lis[[i]]$losing

    gamma.row.lis[[i]] <- .gammaRows_from_comparisons(comp.lis[[i]]$winning,
      losing.mat, sign.vec[i], n)
  }

  mat <- do.call(rbind, gamma.row.lis)
  polyhedra(obj = mat, u = rep(0, nrow(mat)))
}

.gammaRows_from_comparisons <- function(vec, mat, sign.win, n){
  stopifnot(length(vec) == 3, ncol(mat) == 3)

  win.contrast <- .cusum_contrast_full(vec[1], vec[2], vec[3], n)
  lose.contrast <- t(apply(mat, 1, function(x){
    .cusum_contrast_full(x[1], x[2], x[3], n)
  }))

  # add inequalities to compare winning split to all other splits
  res <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win,
    rep(1, nrow(lose.contrast)))
  res2 <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win,
    -rep(1, nrow(lose.contrast)))

  # add inequalities to compare splits to 0 (ensure correct sign)
  rbind(res, res2)
}

.vector_matrix_signedDiff <- function(vec, mat, sign.vec, sign.mat){
  stopifnot(!is.matrix(vec), is.numeric(vec), is.matrix(mat), is.numeric(mat))
  stopifnot(length(vec) == ncol(mat))
  stopifnot(length(sign.vec) == 1, length(sign.mat) == nrow(mat))
  stopifnot(all(c(sign.vec, sign.mat) %in% c(-1,0,1)))

  t(sign.vec * vec - t(sign.mat * mat))
}



##' Function to collect halfsapces for output from the binSeg_fixedThresh(),
##' embedded in an object |obj| of class "bsFt".
##' @param obj object of bsFt class
##' @param verbose Whether or not to print things.
##' @return An object of class polyhedra
polyhedra.bsFt <- function(obj, verbose = F) {

    ## Extracting things
    y = obj$y
    thresh = obj$thresh
    n = length(y)

    ## Initialize G and u
    ii = 1
    G = matrix(NA,nrow = nrow(obj$infotable)*n, ncol = n)
    u = rep(NA, nrow(obj$infotable)*n)
    nrow.G = 0


    ## For each node in Blist, collect halfspaces.
    if(nrow(obj$infotable)==0) return(polyhedra(gamma=G, u=u))
    for(irow in 1:nrow(obj$infotable)){
        myrow = obj$infotable[irow,]
        if(verbose) cat('j,k is', myrow[,"j"], myrow[,"k"], fill=TRUE)


        ## Collect halfspaces
        if(myrow[,"len"]>=1){
            my.halfspaces = halfspaces(s = myrow[,"s"],
                                       b = myrow[,"b"],
                                       e = myrow[,"e"],
                                       z = myrow[,"dir"],
                                       thresh = obj$thresh,
                                       n = n,
                                       y = y,
                                       is.terminal.node = (!myrow[,"pass"]))
        }

        newrows = my.halfspaces[["V"]]
        newconst = my.halfspaces[["u"]]
        ## rownames(newrows)= rep(paste0(j,"-",k),length(newrowinds))

        ## Add to Gamma matrix
        newrowinds = nrow.G + c(1:nrow(newrows))
        if(any(newrowinds > nrow(G) )){
            G = rbind(G, matrix(NA,nrow=nrow(G),ncol=n))
            u = c(u, rep(NA,length(u)))
        }
        G[newrowinds,] = newrows
        u[newrowinds] = newconst

        ## Updates for loop
        nrow.G = nrow.G + nrow(newrows)
        ii = ii + 1
    }
    if(verbose) cat(fill=T)

    ## Trim and return
    G = trim(G,"row")
    u = trim(u)
    return(polyhedra(obj = G, u = u))
}



##' Calculates the halfspace vectors for the maximiiming breakpoint and all the
##' signs, for fixed-threshold SBS.
##' @param is.terminal.node TRUE if no CUSUMs from this node pass the threshold.
##' @param b breakpoint location.
##' @return A list of two objects \code{V} and \code{u}, for the halfspaces in
##'     represented as V'y>u

halfspaces <- function(s, b, e, z, thresh, n, y, is.terminal.node=F , verbose=F){
    if(verbose) cat("s,b,e are", s,b,e, fill=T)
    if(!is.na(b)){
        if(!(s <= b & b <= e)){
            stop("s<=b<=e is not true")
        }
    }

    ## Make empty things, initialize
    V = matrix(NA, nrow = n^2, ncol = n)
    u = rep(NA,n)
    other.bs = (s:(e-1))
    other.bs = other.bs[other.bs != b]
    ii = 0

    ## CUSUM comparison (s,b,e)
    v.this = cusum(s=s, b=b, e=e, y=y,
                   contrast.vec=TRUE, right.to.left=TRUE)
    ## z.this = sign(sum(v.this * y))
    z.this = z
    vz.this = v.this * z.this

    ## 1. Characterizing this break's sign.
    ii = ii+1
    V[ii,] = vz.this
    u[ii] = 0

    ## If it isn't a terminal node
    if(!is.terminal.node){
        ## 2. Characterizing threshold exceedance/nonexceedance
        ii = ii+1
        V[ii,] = vz.this
        u[ii] = thresh

        for(other.b in other.bs){
            v.other = cusum(s=s, b=other.b, e=e, y=y, contrast.vec=T,
                            right.to.left=T)

            ## 3. Characterizing maximizer of |sqrt mean difference|
            if(!is.terminal.node){
                ## other cusum is larger than -|this cusum|
                ii = ii+1
                stopifnot(-sum(vz.this*y) < sum(v.other*y))
                V[ii,] = vz.this + v.other
                u[ii] = 0
                ## other cusum is smaller than +|this cusum|
                ii = ii+1
                stopifnot(sum(v.other*y) < sum(vz.this*y) )
                V[ii,] = vz.this - v.other
                u[ii] = 0
            }
        }
    }

    ## If it is a terminal node
    if(is.terminal.node){
        for(other.b in other.bs){
            v.other = cusum(s=s, b=other.b, e=e, y=y, contrast.vec=T,
                            right.to.left=T)

            ## 4. Characterizing nonexceedance of |cusum|
            if(!is.terminal.node){
                ## other cusum is larger than -thresh
                ii = ii+1
                stopifnot( sum(v.other*y) > -thresh)
                V[ii,] =  v.other + thresh
                u[ii] = thresh
                ## other cusum is smaller than +thresh
                ii = ii+1
                stopifnot(sum(v.other*y) < thresh )
                V[ii,] = - v.other + thresh
                u[ii] = thresh
            }
        }
    }

    V = V[1:ii,,drop=F]
    u = u[1:ii]

    return(list(V=V,u=u))
}

## Testing code
## s = myrow[,"s"];
##                                    b = myrow[,"b"];
##                                    e = myrow[,"e"];
##                                    z = myrow[,"dir"];
##                                    thresh = obj$thresh;
##                                    n = n;
##                                    y = y;
##                                    is.terminal.node = (!myrow[,"pass"])
