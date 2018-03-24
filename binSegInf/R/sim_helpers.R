##' Function to generate one-jump mean
##' @export
onejump <- function(lev, n){
    c(rep(0,n/2),rep(lev,n/2))
}


##' Functions to generate two-jump mean
##' @export
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}


##' Function to generate four-jump mean
##' @export
fourjump <- function(lev,n){c(rep(0,n/5), rep(lev,n/5), rep(0,n/5), rep(-2*lev, n/5), rep(0,n/5) )}
