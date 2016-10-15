##' Constructor for cplist object.
##' @param nrow creates an all-NA Matrix of dimension nrow x 3. Note, it is a
##'     Matrix, not a matrix.
##' @import Matrix
cplist <- function(nrow) {
    structure(list(mat = Matrix(NA, nrow=nrow, ncol=3), last.row=0),
                  class = "cplist")
}

##' Function to add entries to cplist
##' @param cplist list containing an n x 3 matrix, and the index of the last
##'     nonempty (i.e. not all NA's) row.
##' @import Matrix
add.cplist = function(cplist, new.j, new.k, newentry){
    newrow = c(new.j, new.k, newentry)
    cplist$mat[cplist$last.row+1,1:3] = newrow
    cplist$last.row = cplist$last.row+1
    return(cplist)
}

##' Gets the value corresponding to the (j,k)'th entry of Cplist
extract.cplist = function(cplist,j,k){
   return(cplist$mat[where.jk(Cplist,j,k),3])
}

## add <- function(x) UseMethod("add")
## extract <- function(x) UseMethod("extract")

##' Print function
print.cplist = function(cplist){
    if(cplist$last.row==0){ print("Empty cplist object!")
    } else{ print(cplist$mat[1:cplist$last.row,])}
}

## ## Examples
## slist = structure(list(mat = Matrix(NA, nrow=30, ncol=3), last.row=0),
##                   class = "cplist")
## methods(class="cplist")
env = new.env()
env$Blist = cplist(10)
class(env$Blist)
add(env$Blist,2,3,4)
extract(env$Blist,2,3)
env$Blist
