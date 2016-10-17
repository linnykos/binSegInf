##' Constructor for cplist object.
##' @param nrow creates an all-NA Matrix of dimension nrow x 3. Note, it is a
##'     Matrix, not a matrix.
##' @import Matrix
cplist <- function(nrow) {
    structure(list(mat = Matrix(NA, nrow=nrow, ncol=3), last.row=0),
                  class = "cplist")
}


##' Check if object is of class "cplist"
is.cplist <- function(someobj){ inherits(someobj, "cplist") }


add <- function(x) UseMethod("add")
##' Function to add entries to cplist
##' @param cplist list containing an n x 3 matrix, and the index of the last
##'     nonempty (i.e. not all NA's) row.
##' @import Matrix
add.cplist <- function(cplist, new.j, new.k, newentry){

    ## Check if (new.j,new.k) already exists!
    already.exists = !(all(is.na(where_jk.cplist(cplist, new.j, new.k)))) ## This is wrong.
    if(already.exists) stop(paste("j=",new.j, "and", "k=",new.k, "already exist!"))

    ## If cplist$mat is not large enough, then double the size
    numrows = nrow(cplist$mat)
    if(cplist$last.row>= numrows){
        cplist$mat = rbind(cplist$mat, matrix(NA, nrow=numrows, ncol=3))
    }

    ## Append the new row
    newrow = c(new.j, new.k, newentry)
    cplist$mat[cplist$last.row+1,1:3] = newrow
    cplist$last.row = cplist$last.row+1
    return(cplist)
}

extract <- function(x) UseMethod("extract")
#' Gets the value corresponding to the (j,k)'th entry of Cplist
extract.cplist <- function(cplist,j,k){
   return(cplist$mat[where_jk.cplist(cplist,j,k),3])
}

##' Print function
print.cplist <- function(cplist){
    ## if(cplist$last.row==0){ print("Empty cplist object!")
    ## } else{ print(cplist$mat[1:cplist$last.row,])}
    ## print(trim.cplist(cplist))
    print(cplist$mat)
    ## print(cplist$last.row)
}


trim <- function(x) UseMethod("trim")
##' Trim function
trim.cplist <- function(cplist){
    if(cplist$last.row==0){
        return(rbind(c(1,1,1))[-1,,drop=FALSE])
    } else {
       return(cplist$mat[1:cplist$last.row,,drop=FALSE])
    }
}



where_jk <- function(x) UseMethod("where_jk")
##' Search in cplist$mat for the couplet (j,k) in the first two columns e.g. if
##' j=13 and k = 39, then it searches for the row (13,39,XXX) in an n by 3
##' matrix.
where_jk.cplist <- function(cplist, j, k, warn=FALSE){
    ## cat("j,k queried for are", j, k, fill=TRUE)
    
    ## Basic check
    stopifnot(is.cplist(cplist))

    # Extract columns containing j and k.
    jvec = trim.cplist(cplist)[,1]
    kvec = trim.cplist(cplist)[,2]

    ## ## Find any thing that matches j
    ## one.j.loc = binary_search(jvec, goal=j)
    ## ## From there, crawl to find /all/ that matches j
    ## all.j.loc = crawler(jvec, one.j.loc, j)

    all.j.loc = which(jvec==j)
    if(length(all.j.loc)==0){
        if(warn) warning(paste("The indices (",j,",",k,")", "don't exist!"))
        return(NA)
    }

    ## From there, crawl to find /all/ that matches k
    ## k.loc = min(all.j.loc)-1 + binary_search(kvec[all.j.loc], goal=k)
    k.loc = min(all.j.loc)-1 + which(kvec[all.j.loc]==k)
    if(length(k.loc)==0){#kvec[k.loc] != k){
        if(warn) warning(paste("The indices (",j,",",k,")", "don't exist!"))
        return(NA)
    }

    ## Return the location.
    return(k.loc)
}



## methods(class="cplist")
## a=cplist(4)
## class(a)
## a=add.cplist(a,1,1,3)
## a=add.cplist(a,2,3,3)
## a=add.cplist(a,3,5,3)
## where_jk.cplist(a,4,4)
## where_jk(a,4,4)
## add(cplist(1),1,1,1)

    ## stopifnot(all(jvec[all.j.loc]==j) & all(jvec[!all.j.loc]!=j) )
