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



##' Function to collect polyhedra given some output from the
##' binSeg_fixedThresh(), of class "bsFt".
##' 
##' @param obj object of bsFt class
##' @param verbose Whether or not to print things.
##'
##' @return An object of class polyhedra
polyhedra.bsFt = function(obj, verbose = F) {

    ## Extracting things
    Blist = obj$bs.output$Blist
    blist = obj$bs.output$blist
    slist = obj$bs.output$slist
    elist = obj$bs.output$elist
    selist = obj$bs.output$selist
    zlist = obj$bs.output$zlist
    Zlist = obj$bs.output$Zlist
    y = obj$y
    thresh = obj$thresh
    n = length(y)

    ## Initialize G and u
    ii = 1
    G = matrix(NA,nrow = nrow((Blist))*n, ncol = n)
    u = rep(NA, nrow((Blist))*n)
    nrow.G = 0


    ## For each node in Blist, collect halfspaces.
    for(my.j in 1:nrow((Blist))){
        b = Blist[my.j, "val"]
        z = Zlist[my.j, "val"]
        j = Blist[my.j,"j"]-1
        k = Blist[my.j,"k"]
        if(verbose) cat('j,k is', j,k, fill=TRUE)
        ## Extract s,e, check if terminal node
        ind = which(selist[,"j"]==j&selist[,"k"]==k)
        se = c(slist[ind,"val"], elist[ind,"val"])
        is.terminal.node = !any(apply(blist[,c(1,2)],1,function(myrow)(myrow[1]==j+1)&&(myrow[2]==k)))

        ## Collect halfspaces
        my.halfspaces = halfspaces(s = se[1],
                                   e = se[2],
                                   b = b,
                                   z=z,
                                   thresh = thresh,
                                   n = n,
                                   y = y,
                                   is.terminal.node = is.terminal.node)

        newrows = my.halfspaces[["V"]]
        newconst = my.halfspaces[["u"]]
        
        ## Move on if no comparisons were made i.e. lengths to left=2, right=1
        if(dim(newrows)[1]==0) next
            
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
