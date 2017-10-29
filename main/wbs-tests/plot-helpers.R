##' Add title according to simulation setting
maketitle <- function(lev=NULL, n=NULL, other=NULL){
    assert_that(!(is.null(lev)&is.null(n)))
    mytitle = ""
    if(!is.null(lev)) mytitle = paste0(mytitle, "Jump size = ", lev)
    if(!is.null(n)) mytitle = paste0(mytitle, " n = ", n)
    if(!is.null(other)) mytitle = paste(mytitle, "\n", othertitle)
    title(main=mytitle)
}
