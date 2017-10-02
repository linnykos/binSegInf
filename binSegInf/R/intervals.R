##' Constructor for intervals object.
##' @param nrow creates an all-NA matrix of dimension nrow x 3. The first two
##'     columns must be the numeric (no check yet), but the last column can be
##'     of any type you want. Initializes to numeric.
##' @export
intervals <- function(nrow) {

    ## generate_intervals(n, numIntervals, seed=NULL, start.end.list = NULL){
    structure(list(starts=starts,
                   ends=ends,
                   intervals=intervals,
                   se=Map(c,starts,ends)),
              class="intervals")
}

##' Checks if interval is a valid object of the |interval| class.
is_valid.intervals <- function(obj){
    return(all(names(intervals) %in% c("starts", "ends", "intervals", "se")))
}



##' After making intervals, you can attempt to plot them.
##' @export
plot.intervals <- function(intervals){
    ## Basic checks
    stopifnot(is_valid.intervals(intervals))

    numIntervals = length(intervals$se)
    graphics::plot(NA,
                ylim = c(0,numIntervals),
                xlim = c(0,max(intervals$e)),
                xlab = "intervals",
                ylab = "")
    for(ii in 1:numIntervals){
        se = intervals$se[[ii]]
        graphics::lines(x=se, y = c(ii,ii))
    }
}


##' Add a single interval to existing set of intervals
##' @param old.intervals an object of class |intervals|
##' @param new.s new start point
##' @param new.e new end point
##' @export
add.intervals <- function(intervals, new.s, new.e){
    ## Basic checks
    stopifnot(new.s < new.e)

    ## Append new interval information
    new.starts = c(intervals$starts, new.s)
    new.ends = c(intervals$ends, new.e)
    new.intervals = c(intervals$intervals, list(new.s:new.e))

    ## Return as |intervals| class object
    structure(list(starts=new.starts,
                   ends=new.ends,
                   intervals=new.intervals,
                   se=Map(c,new.starts,new.ends)),
              class="intervals")
}
