## name = "onejump"
name = "dummy"
## setwd(paste0("~/Dropbox/teaching/f15-350/lectures/",name,"/"))
rmarkdown::render(paste0(name,".Rmd"),output_format="pdf_document",output_file=paste0(name,".pdf"),quiet=TRUE)
## rmarkdown::render(paste0(name,".Rmd"),output_format="slidy_presentation",output_file=paste0(name,"-pres.html"),quiet=TRUE)
