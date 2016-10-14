  #filename = "./onejump.Rmd"
  #filename = "./twojump.Rmd"
  filename = "./threejump.Rmd"
  library(rmarkdown)
  render(filename, html_document(toc=TRUE))