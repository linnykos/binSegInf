##' Add title according to simulation setting
maketitle <- function(lev=NULL, n=NULL, othertitle=NULL){
    assert_that(!(is.null(lev)&is.null(n)))
    mytitle = ""
    if(!is.null(lev)) mytitle = paste0(mytitle, "Jump size = ", lev)
    if(!is.null(n)) mytitle = paste0(mytitle, " n = ", n)
    if(!is.null(othertitle)) mytitle = paste(mytitle, "\n", othertitle)
    title(main=mytitle)
}




plot.visc <- function(pvs,locs,visc,n,levs, mytitle="One jump, one step",
                      prefix = "jump-size="){
    cond.pvs <- Map(function(mypv, myloc){
        mypv[myloc %in% visc]},pvs, locs)
    names(cond.pvs) = paste(prefix, levs)
    ## rid.of.ones=TRUE
    ## if(rid.of.ones) cond.pvs = lapply(cond.pvs, function(mypvs)mypvs[which(mypvs!=1)])
    qqunif(cond.pvs,cols=1:4, xlab="", ylab="")
    maketitle(n=n, othertitle=mytitle)
}

plot.nulls <- function(pvs,locs,truths,n,levs, mytitle="One jump, three steps"){
    cond.pvs = Map(function(mytruth, mypv)mypv[mytruth], truths, pvs)
    names(cond.pvs) = paste("jump-size=",levs)
    ## rid.of.ones=TRUE
    ## if(rid.of.ones) cond.pvs = lapply(cond.pvs, function(mypvs)mypvs[which(mypvs!=1)])
    qqunif(cond.pvs,cols=1:4, xlab="", ylab="")
    maketitle(n=n, othertitle=mytitle)
}
